# News

## Version 0.1.0

- Fixed registered-package backend builds by moving Julia-side build discovery
  into `deps/build.jl`; the CMake build no longer runs Julia with the package
  source directory as its project.
- Moved the generated shared library to a package scratch space so locally
  compiled backends survive package source updates.
- Stopped exporting `log` so `using GridapTrilinos` does not shadow
  `Base.log`; use the `solver.log` property for solver results.
- Declared the package test target so `Pkg.test("GridapTrilinos")` includes
  test-only dependencies, and switched CI to `julia-actions/julia-runtest`.
- Added package version metadata for Julia package registration.
- Made installation succeed without a local Trilinos installation. The package
  can be loaded without the optional C++ backend, and Trilinos-dependent calls
  report how to build it.
- Moved the generated shared library from `src/GridapTrilinos.so` to
  `deps/usr/lib/GridapTrilinos.so`.
- Clarified the public API as `TrilinosSolve`, `SolverResult`, and the
  solver-result accessors.
- Added Documenter.jl documentation and docstrings for the public API.
- Added a GitHub Actions documentation workflow.
- Added TagBot release automation.
- Added Codecov/code-coverage workflow configuration.
- Added GitHub repository description and topics for package discovery.
