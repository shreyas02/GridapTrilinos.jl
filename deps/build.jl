using Libdl
using MPIPreferences
using Scratch
using libcxxwrap_julia_jll

const GRIDAPTRILINOS_UUID = Base.UUID("fd8928cd-657a-5e46-9bb5-36bf17cf2d1c")
const SCRATCH_KEY = "trilinos-backend"

source_dir = normpath(@__DIR__, "..", "src", "Sharedlib")

if !isdir(source_dir)
  error("GridapTrilinos C++ source directory was not found at $(source_dir)")
end

if !haskey(ENV, "TRILINOS_ROOT") || isempty(ENV["TRILINOS_ROOT"])
  @warn(
    "TRILINOS_ROOT is not set. Set it to your Trilinos installation " *
    "prefix, then run `Pkg.build(\"GridapTrilinos\")` to enable " *
    "Trilinos solves. Skipping the optional C++ build.",
  )
  exit()
end

trilinos_root = ENV["TRILINOS_ROOT"]

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

scratch_root = get_scratch!(GRIDAPTRILINOS_UUID, SCRATCH_KEY, GRIDAPTRILINOS_UUID)
build_dir = joinpath(scratch_root, "build")
install_dir = joinpath(scratch_root, "usr")
lib_dir = joinpath(install_dir, "lib")

mkpath(build_dir)
mkpath(lib_dir)

trilinos_lib_paths = String[]
for lib_subdir in ("lib", "lib64")
  path = joinpath(trilinos_root, lib_subdir)
  isdir(path) && push!(trilinos_lib_paths, path)
end

if !isempty(trilinos_lib_paths)
  ENV["LD_LIBRARY_PATH"] = join(
    vcat(trilinos_lib_paths, get(ENV, "LD_LIBRARY_PATH", "") == "" ? String[] : [ENV["LD_LIBRARY_PATH"]]),
    ":",
  )
end

ENV["CMAKE_PREFIX_PATH"] = join(
  vcat([trilinos_root], get(ENV, "CMAKE_PREFIX_PATH", "") == "" ? String[] : [ENV["CMAKE_PREFIX_PATH"]]),
  ":",
)

trilinos_dir_lib = joinpath(trilinos_root, "lib", "cmake", "Trilinos")
trilinos_dir_lib64 = joinpath(trilinos_root, "lib64", "cmake", "Trilinos")
if isdir(trilinos_dir_lib)
  ENV["Trilinos_DIR"] = trilinos_dir_lib
elseif isdir(trilinos_dir_lib64)
  ENV["Trilinos_DIR"] = trilinos_dir_lib64
end

julia_include_dir = normpath(Sys.BINDIR, "..", "include", "julia")
julia_libjulia = Libdl.dlpath("libjulia")

cmake_configure = `cmake
  -S $(source_dir)
  -B $(build_dir)
  -DGRIDAPTRILINOS_LIBRARY_OUTPUT_DIR=$(lib_dir)
  -DLIBCXXWRAP_ARTIFACT_DIR=$(libcxxwrap_julia_jll.artifact_dir)
  -DJULIA_INCLUDE_DIR=$(julia_include_dir)
  -DJULIA_LIBJULIA=$(julia_libjulia)`

run(cmake_configure)
run(`cmake --build $(build_dir)`)
