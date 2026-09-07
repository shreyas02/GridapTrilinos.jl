# GridapTrilinos.jl

GridapTrilinos provides a `Gridap.Algebra.LinearSolver` implementation backed by
Trilinos. The Julia package can be installed and loaded without Trilinos; calls
that need the optional C++ backend report how to build it.

## Public API

```@docs
GridapTrilinos.TrilinosSolve
GridapTrilinos.SolverResult
GridapTrilinos.name
GridapTrilinos.num_iters
GridapTrilinos.residual
GridapTrilinos.solve_time
GridapTrilinos.verbose
GridapTrilinos.depth
```

## Backend

Build the optional backend after configuring MPI.jl to use the same system MPI
as Trilinos. For a registered package installation:

```bash
export TRILINOS_ROOT=/path/to/TrilinosInstall
julia -e 'using Pkg; Pkg.build("GridapTrilinos")'
```

For a source checkout with the repository environment activated, use:

```bash
export TRILINOS_ROOT=/path/to/TrilinosInstall
julia --project=. -e 'using Pkg; Pkg.build("GridapTrilinos")'
```

The shared library is generated in Julia scratch space:

```text
DEPOT_PATH[1]/scratchspaces/<GridapTrilinos uuid>/trilinos-backend/usr/lib/GridapTrilinos.so
```

Because the backend is outside the versioned package source directory, it is not
lost when Julia replaces the package source tree during an update.
