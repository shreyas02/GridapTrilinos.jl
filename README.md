# GridapTrilinos

GridapTrilinos is a Julia interface for using Trilinos linear solvers from
Gridap/GridapDistributed workflows. The Julia side provides a `TrilinosSolve`
linear solver, while the Trilinos, Tpetra, Belos, FROSch, and Kokkos calls live
in a small C++ shared library exposed to Julia with CxxWrap.

It is authored by Shreyas Prashanth.

## Workflow

There are four parts to the workflow:

1. Install the Julia package dependencies.
2. Configure `MPI.jl` to use the same MPI implementation as Trilinos.
3. Build the C++ shared library against your local Trilinos installation.
4. Load `GridapTrilinos` and use `TrilinosSolve` as a Gridap linear solver.

## Julia Setup

From the repository root:

```julia
using Pkg
Pkg.activate(".")
Pkg.instantiate()
```

This installs the Julia dependencies, including `Gridap`, `GridapDistributed`,
`PartitionedArrays`, `CxxWrap`, `MPIPreferences`, and `libcxxwrap_julia_jll`.

## MPI Compatibility

The MPI used by Julia must be the same MPI implementation used to build
Trilinos. This package uses Trilinos' MPI on the C++ side, while
`GridapDistributed` reaches MPI through `MPI.jl`. If those are different MPI
runtimes, the process can fail at load time or hang/crash during communication.

Configure `MPI.jl` through `MPIPreferences` before building this package. For a
typical system MPI installation:

```julia
using Pkg
Pkg.activate(".")
Pkg.add("MPIPreferences")

using MPIPreferences
MPIPreferences.use_system_binary()
```

If `libmpi` is not in the dynamic linker search path, pass the library location:

```julia
using MPIPreferences
MPIPreferences.use_system_binary(;
    library_names = ["/path/to/mpi/lib/libmpi.so"],
    mpiexec = "/path/to/mpi/bin/mpiexec",
)
```

Restart Julia after changing `MPIPreferences`. Then instantiate/build the
package. `deps/build.jl` checks that `MPI.jl` is configured for `"system"` MPI
before compiling the C++ wrapper.

## Building The C++ Library

The C++ wrapper must be built before calling the Trilinos solver. The generated
file is a shared library, not a portable source file, so each user builds it
locally against their own Julia, MPI, and Trilinos installation.

Set `TRILINOS_ROOT` to the Trilinos installation prefix, then build through
Julia's package build step:

```bash
export TRILINOS_ROOT=/path/to/TrilinosInstall
julia --project=. -e 'using Pkg; Pkg.build("GridapTrilinos")'
```

This calls `deps/build.jl`, which delegates to `src/Sharedlib/configure.sh`.
You can also call the script directly:

```bash
export TRILINOS_ROOT=/path/to/TrilinosInstall
src/Sharedlib/configure.sh
```

The build configures CMake in `src/Sharedlib/build/` and creates:

```text
src/GridapTrilinos.so
```

This file is generated and ignored by git. If it is missing, `using
GridapTrilinos` still works, but any call into the Trilinos solver throws an
error telling you to build the shared library first.

### Custom Trilinos Solve Source

By default, the shared library builds the solver implementation in:

```text
src/Sharedlib/TrilinosSolve.cpp
```

Advanced users can provide a different C++ source file at build time with
`GRIDAPTRILINOS_SOLVE_SOURCE`:

```bash
export TRILINOS_ROOT=/path/to/TrilinosInstall
export GRIDAPTRILINOS_SOLVE_SOURCE=/path/to/MyTrilinosSolve.cpp
src/Sharedlib/configure.sh
```

Relative paths are resolved from `src/Sharedlib/`, so this also works:

```bash
export TRILINOS_ROOT=/path/to/TrilinosInstall
export GRIDAPTRILINOS_SOLVE_SOURCE=MyTrilinosSolve.cpp
src/Sharedlib/configure.sh
```

The custom source must define the same `TrilinosSolve` function declared in
`src/Sharedlib/TrilinosTypes.hpp`:

```cpp
TrilinosSolveData TrilinosSolve(
  const RCP<crs_matrix_type>& A,
  const RCP<vec_type>& b,
  const std::string& parameterFilePath,
  bool verbose);
```

The Gridap/Tpetra construction and Julia wrapper code still comes from
`src/Sharedlib/TrilinosInterface.cpp`; only the Trilinos solve implementation is
replaced.

## Initialisation

Loading the package performs the CxxWrap and Kokkos setup:

```julia
using GridapTrilinos
```

When `src/GridapTrilinos.so` exists, the package:

- loads the C++ module with CxxWrap,
- calls `KokkosInitialize()` during Julia module initialisation,
- registers `KokkosFinalize()` with `atexit`.

You normally do not need to call `KokkosInitialize()` or `KokkosFinalize()`
manually.

## Usage

Create a solver from a Trilinos XML parameter file:

```julia
using GridapTrilinos

solver = TrilinosSolve("path/to/trilinos_parameters.xml")
```

Optionally set the estimated maximum number of nonzeros per local matrix row:

```julia
solver = TrilinosSolve("path/to/trilinos_parameters.xml", 100)
```

`TrilinosSolve` implements the Gridap linear solver interface, so it is intended
to be passed wherever a `Gridap.Algebra.LinearSolver` is expected. Internally,
numerical setup builds and stores the distributed Tpetra matrix. During
`solve!`, the package builds the Tpetra right-hand side, calls the wrapped
Trilinos solve implementation, and copies the local solution back into the
Gridap vector.

After a solve has run, inspect the recorded solver result:

```julia
result = solver.log

result.name
result.num_iters
result.residual
result.solve_time
```

The lower-level C++ wrapper functions are implementation details; typical
Gridap usage should go through `TrilinosSolve`.

## Development Checks

Run the package tests from the repository root:

```bash
julia --project=. test/runtests.jl
```

The default test run checks that the distributed Poisson workflow is loadable,
but skips the full MPI/Trilinos solve. To run the solver workflow test, add the
XML file at `test/poisson_frosch.xml` and launch with 4 MPI ranks:

```bash
GRIDAPTRILINOS_RUN_SOLVER_TEST=true \
  mpiexecjl --project=. -n 4 julia test/runtests.jl
```

To use a different XML file:

```bash
GRIDAPTRILINOS_RUN_SOLVER_TEST=true \
GRIDAPTRILINOS_PARAMETER_FILE=test/my_parameters.xml \
  mpiexecjl --project=. -n 4 julia test/runtests.jl
```

Rebuild the C++ library after changing files in `src/Sharedlib/`:

```bash
src/Sharedlib/configure.sh
```
