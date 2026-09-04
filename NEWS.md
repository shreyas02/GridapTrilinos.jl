# News

## Version 0.1.0

- Added package version metadata for Julia package registration.
- Made installation succeed without a local Trilinos installation. The package
  can be loaded without the optional C++ backend, and Trilinos-dependent calls
  report how to build it.
- Moved the generated shared library from `src/GridapTrilinos.so` to
  `deps/usr/lib/GridapTrilinos.so`.
- Clarified the public API as `TrilinosSolve`, `SolverResult`, `log`, and the
  solver-result accessors.
- Added Documenter.jl documentation and docstrings for the public API.
- Added a GitHub Actions documentation workflow.
- Added TagBot release automation.
- Added Codecov/code-coverage workflow configuration.
- Added GitHub repository description and topics for package discovery.
