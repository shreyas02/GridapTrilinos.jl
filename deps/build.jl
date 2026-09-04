build_script = normpath(@__DIR__, "..", "src", "Sharedlib", "configure.sh")

using MPIPreferences

if !isfile(build_script)
    error("GridapTrilinos build script was not found at $(build_script)")
end

if !haskey(ENV, "TRILINOS_ROOT") || isempty(ENV["TRILINOS_ROOT"])
  @warn(
    "TRILINOS_ROOT is not set. Set it to your Trilinos installation " *
    "prefix, then run `Pkg.build(\"GridapTrilinos\")` to enable " *
    "Trilinos solves. Skipping the optional C++ build.",
  )
  exit()
end

if MPIPreferences.binary != "system"
  @warn(
    "MPI.jl is configured to use $(MPIPreferences.binary), but " *
    "GridapTrilinos must use the same MPI as the Trilinos installation. " *
    "Configure MPI.jl to use the system MPI from your Trilinos/MPI stack " *
    "with `MPIPreferences.use_system_binary(...)`, restart Julia, then " *
    "run `Pkg.build(\"GridapTrilinos\")` again. Skipping the optional C++ " *
    "build.",
  )
  exit()
end

run(`$(build_script)`)
