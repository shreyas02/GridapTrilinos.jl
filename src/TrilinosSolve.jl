#######################
## Begin Trilinos Solve 
struct TrilinosSolve <: Gridap.Algebra.LinearSolver
    parameter_file::String
    log::Base.RefValue{Union{Nothing,SolverResultAllocated}}
end

function TrilinosSolve(parameter_file::AbstractString)
    return TrilinosSolve(
        String(parameter_file),
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

mutable struct TrilinosSolveNumericalSetup{M,R} <: Gridap.Algebra.NumericalSetup
    solver::TrilinosSolve
    matrix_partition::M
    row_partition::R
end

function Gridap.Algebra.symbolic_setup(solver::TrilinosSolve, mat::AbstractMatrix)
    return TrilinosSolveSymbolicSetup(solver)
end

function Gridap.Algebra.numerical_setup(ss::TrilinosSolveSymbolicSetup, A::AbstractMatrix)
    matrix_partition = _construct_tpetra_matrix_partition(A)
    return TrilinosSolveNumericalSetup(ss.solver, matrix_partition, A.row_partition)
end

function Gridap.Algebra.numerical_setup(
    ss::TrilinosSolveSymbolicSetup,
    A::AbstractMatrix,
    x::AbstractVector,
)
    matrix_partition = _construct_tpetra_matrix_partition(A)
    return TrilinosSolveNumericalSetup(ss.solver, matrix_partition, A.row_partition)
end

function Gridap.Algebra.numerical_setup!(ns::TrilinosSolveNumericalSetup, A::AbstractMatrix)
    matrix_partition = _construct_tpetra_matrix_partition(A)
    ns.matrix_partition = matrix_partition
    ns.row_partition = A.row_partition
    return ns
end

function _construct_tpetra_matrix_partition(A::AbstractMatrix)
    map(A.matrix_partition, A.row_partition, A.col_partition) do LocMat, RowMap, ColMap
        OwnToLocalRow = own_to_local(RowMap)
        ConstructTpetraMatrixWrapper(
            LocMat.nzval,
            LocMat.colval,
            LocMat.rowptr,
            collect(local_to_global(RowMap)),
            collect(local_to_global(ColMap)),
            size(A)[1],
            length(OwnToLocalRow),
            OwnToLocalRow,
        )
    end
end

function Gridap.Algebra.solve!(
    x::AbstractVector,
    ns::TrilinosSolveNumericalSetup,
    b::AbstractVector,
)
    map(
        ns.matrix_partition,
        ns.row_partition,
        b.vector_partition,
        x.vector_partition,
        x.index_partition,
    ) do TpetraMatrix, RowMap, LocRhs, LocSoln, xmap
        ## # Construct Tpetra RHS vector
        OwnToLocalRow = own_to_local(RowMap)
        OwnToLocalSol = own_to_local(xmap)
        TpetraRhs = ConstructTpetraVectorWrapper(
            LocRhs,
            length(OwnToLocalRow),
            OwnToLocalRow,
            TpetraMatrix,
        )
        ## Solve the linear system using Trilinos
        solve_data = TrilinosSolveWrapper(
            TpetraMatrix,
            TpetraRhs,
            ns.solver.parameter_file,
        )
        ## Copy the solution back to the Gridap vector
        getfield(ns.solver, :log)[] = CopySolutionWrapper(
            solve_data,
            LocSoln,
            length(OwnToLocalRow),
            OwnToLocalSol,
        )
    end
    return ns
end
## End Trilinos Solve
#####################
