#######################
## Begin Trilinos Solve 
struct TrilinosSolve <: Gridap.Algebra.LinearSolver
    parameter_file::String
    max_entries_per_row::Int
    log::Base.RefValue{Union{Nothing,SolverResultAllocated}}
end

function TrilinosSolve(parameter_file::AbstractString, max_entries_per_row::Integer=100)
    return TrilinosSolve(
        String(parameter_file),
        Int(max_entries_per_row),
        Ref{Union{Nothing,SolverResultAllocated}}(nothing),
    )
end

log(solver::TrilinosSolve) = getfield(solver, :log)[]

function Base.getproperty(solver::TrilinosSolve, name::Symbol)
    if name === :log
        value = getfield(solver, :log)[]
        value === nothing && error("Solver log is not available yet. Run the solve first.")
        return value
    end
    return getfield(solver, name)
end

struct TrilinosSolveSymbolicSetup <: Gridap.Algebra.SymbolicSetup
    solver::TrilinosSolve
end

mutable struct TrilinosSolveNumericalSetup{T<:AbstractMatrix} <: Gridap.Algebra.NumericalSetup
    solver::TrilinosSolve
    A::T
end

function Gridap.Algebra.symbolic_setup(solver::TrilinosSolve, mat::AbstractMatrix)
    return TrilinosSolveSymbolicSetup(solver)
end

function Gridap.Algebra.numerical_setup(ss::TrilinosSolveSymbolicSetup, A::AbstractMatrix)
    return TrilinosSolveNumericalSetup(ss.solver, A)
end

function Gridap.Algebra.numerical_setup(
    ss::TrilinosSolveSymbolicSetup,
    A::AbstractMatrix,
    x::AbstractVector,
)
    return TrilinosSolveNumericalSetup(ss.solver, A)
end

function Gridap.Algebra.numerical_setup!(ns::TrilinosSolveNumericalSetup, A::AbstractMatrix)
    return ns
end

function Gridap.Algebra.solve!(
    x::AbstractVector,
    ns::TrilinosSolveNumericalSetup,
    b::AbstractVector,
)
    map(
        ns.A.matrix_partition,
        ns.A.row_partition,
        ns.A.col_partition,
        b.vector_partition,
        x.vector_partition,
        x.index_partition,
    ) do LocMat, RowMap, ColMap, LocRhs, LocSoln, xmap
        OwnRowMap = RowMap[own_to_local(RowMap)]
        getfield(ns.solver, :log)[] = TrilinosParallel(
            LocMat.nzval,
            LocMat.rowval .- 1,
            LocMat.colptr .- 1,
            LocRhs,
            LocSoln,
            OwnRowMap .- 1,
            ColMap .- 1,
            size(ns.A)[1],
            size(OwnRowMap)[1],
            size(LocMat)[2],
            own_to_local(xmap) .- 1,
            own_to_local(RowMap) .- 1,
            ns.solver.max_entries_per_row,
            ns.solver.parameter_file,
        )
    end
    return ns
end
## End Trilinos Solve
#####################
