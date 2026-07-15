# Load the module and generate the functions
module GridapTrilinos
  using Gridap
  using CxxWrap
  using Libdl
  using MPIPreferences
  import PartitionedArrays: local_to_global, own_to_local

  # cpp structs and functions to be exported to julia
  export KokkosInitialize, KokkosFinalize
  export SolverResult, SolverResultAllocated, SolverResultDereferenced
  export TrilinosSolve
  export log
  # properties of the SolverResult struct to be exported to julia
  export num_iters, residual, solve_time, name, verbose, depth

  const _sharedlib_path = joinpath(@__DIR__, "GridapTrilinos")
  const _sharedlib_file = _sharedlib_path * "." * Libdl.dlext
  const _has_sharedlib = isfile(_sharedlib_file) && MPIPreferences.binary == "system"

  if _has_sharedlib
    @wrapmodule(() -> _sharedlib_path)
  else
    struct SolverResult end
    struct SolverResultAllocated end
    struct SolverResultDereferenced end

    function _missing_sharedlib()
      if isfile(_sharedlib_file) && MPIPreferences.binary != "system"
        error(
          "GridapTrilinos found $(_sharedlib_file), but MPI.jl is configured to use " *
          "$(MPIPreferences.binary). Configure MPI.jl to use the same system MPI as " *
          "Trilinos with `MPIPreferences.use_system_binary(...)`, restart Julia, " *
          "then rebuild GridapTrilinos.",
        )
      end
      error(
        "GridapTrilinos shared library was not found at $(_sharedlib_file). " *
        "Build it first with `src/Sharedlib/configure.sh` after setting TRILINOS_ROOT.",
      )
    end

    ConstructTpetraMatrixWrapper(args...) = _missing_sharedlib()
    ConstructTpetraVectorWrapper(args...) = _missing_sharedlib()
    TrilinosSolveWrapper(args...) = _missing_sharedlib()
    CopySolutionWrapper(args...) = _missing_sharedlib()
    KokkosInitialize() = nothing
    KokkosFinalize() = nothing
    num_iters(result::Union{SolverResultAllocated,SolverResultDereferenced}) = _missing_sharedlib()
    residual(result::Union{SolverResultAllocated,SolverResultDereferenced}) = _missing_sharedlib()
    solve_time(result::Union{SolverResultAllocated,SolverResultDereferenced}) = _missing_sharedlib()
    name(result::Union{SolverResultAllocated,SolverResultDereferenced}) = _missing_sharedlib()
    verbose(result::Union{SolverResultAllocated,SolverResultDereferenced}) = _missing_sharedlib()
    depth(result::Union{SolverResultAllocated,SolverResultDereferenced}) = _missing_sharedlib()
  end

  const WrappedSolverResult = Union{SolverResultAllocated,SolverResultDereferenced}

  const _solver_result_properties = (:num_iters, :residual, :solve_time, :name, :verbose, :depth)

  function Base.getproperty(result::WrappedSolverResult, element::Symbol)
    if element === :num_iters
      return num_iters(result)
    elseif element === :residual
      return residual(result)
    elseif element === :solve_time
      return solve_time(result)
    elseif element === :name
      return name(result)
    elseif element === :verbose
      return verbose(result)
    elseif element === :depth
      return depth(result)
    end
    return getfield(result, element)
  end

  Base.propertynames(::WrappedSolverResult, private::Bool=false) =
    private ? (_solver_result_properties..., fieldnames(SolverResultAllocated)...) : _solver_result_properties

  include("TrilinosSolve.jl")

  function __init__()
    if _has_sharedlib
      @initcxx
      KokkosInitialize()
      atexit(KokkosFinalize)
    end
  end
end
