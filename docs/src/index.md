# GridapTrilinos.jl

GridapTrilinos provides a `Gridap.Algebra.LinearSolver` implementation backed by
Trilinos. The Julia package can be installed and loaded without Trilinos; calls
that need the optional C++ backend report how to build it.

## Public API

```@docs
GridapTrilinos.TrilinosSolve
GridapTrilinos.SolverResult
GridapTrilinos.log
GridapTrilinos.name
GridapTrilinos.num_iters
GridapTrilinos.residual
GridapTrilinos.solve_time
GridapTrilinos.verbose
GridapTrilinos.depth
```

## Backend

Build the optional backend after configuring MPI.jl to use the same system MPI
as Trilinos:

```bash
export TRILINOS_ROOT=/path/to/TrilinosInstall
julia --project=. -e 'using Pkg; Pkg.build("GridapTrilinos")'
```

The shared library is generated under:

```text
deps/usr/lib/GridapTrilinos.so
```
