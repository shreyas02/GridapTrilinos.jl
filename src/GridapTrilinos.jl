# Load the module and generate the functions
module GridapTrilinos
  using Gridap
  using CxxWrap
  using Libdl
  using MPIPreferences
  import PartitionedArrays: local_to_global, own_to_local

  export TrilinosSolve, SolverResult, log
  export num_iters, residual, solve_time, name, verbose, depth

  const _sharedlib_path = joinpath(dirname(@__DIR__), "deps", "usr", "lib", "GridapTrilinos")
  const _sharedlib_file = _sharedlib_path * "." * Libdl.dlext

  function _probe_sharedlib()
    isfile(_sharedlib_file) || return (false, nothing)
    MPIPreferences.binary == "system" || return (false, nothing)
    handle = Libdl.dlopen_e(_sharedlib_file)
    if handle == C_NULL
      return (false, Libdl.dlerror())
    end
    Libdl.dlclose(handle)
    return (true, nothing)
  end

  const _sharedlib_probe = _probe_sharedlib()
  const _has_sharedlib = first(_sharedlib_probe)
  const _sharedlib_load_error = last(_sharedlib_probe)

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
      if isfile(_sharedlib_file) && _sharedlib_load_error !== nothing
        error(
          "GridapTrilinos found $(_sharedlib_file), but it could not be loaded: " *
          "$(_sharedlib_load_error). Ensure the Trilinos and MPI shared libraries " *
          "used to build GridapTrilinos are available in the dynamic linker path.",
        )
      end
      error(
        "GridapTrilinos shared library was not found at $(_sharedlib_file). " *
        "Set TRILINOS_ROOT and run `Pkg.build(\"GridapTrilinos\")` to enable Trilinos solves.",
      )
    end

    ConstructTpetraMatrixWrapper(args...) = _missing_sharedlib()
    ConstructTpetraVectorWrapper(args...) = _missing_sharedlib()
    TrilinosSolverSetupWrapper(args...) = _missing_sharedlib()
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

  """
      SolverResult

  Result object returned by the Trilinos backend after a solve has completed.

  Access values with the public properties `name`, `num_iters`, `residual`,
  `solve_time`, `verbose`, and `depth`, or with the matching accessor
  functions. A solver's latest result is available as `solver.log` after
  `Gridap.Algebra.solve!` has run.
  """
  SolverResult

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

  """
      log(solver::TrilinosSolve)

  Return the latest `SolverResult` recorded by `solver`, or `nothing` if no
  solve has completed.

  The property form `solver.log` is stricter and throws when no result is
  available.
  """
  log

  """
      name(result::SolverResult)

  Return the Trilinos linear solver and preconditioner name reported for a
  completed solve.
  """
  name

  """
      num_iters(result::SolverResult)

  Return the iteration count reported by Trilinos for a completed solve.
  """
  num_iters

  """
      residual(result::SolverResult)

  Return the achieved residual tolerance reported by Trilinos for a completed
  solve.
  """
  residual

  """
      solve_time(result::SolverResult)

  Return the elapsed Trilinos solve time in seconds.
  """
  solve_time

  """
      verbose(result::SolverResult)

  Return the backend verbosity value stored in a `SolverResult`.
  """
  verbose

  """
      depth(result::SolverResult)

  Return the backend depth value stored in a `SolverResult`.
  """
  depth

  function __init__()
    if _has_sharedlib
      @initcxx
      KokkosInitialize()
      atexit(KokkosFinalize)
    end
  end
end
