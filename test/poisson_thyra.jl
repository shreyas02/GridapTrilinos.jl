using Gridap
using GridapDistributed
using GridapTrilinos
using MPI
using PartitionedArrays
using SparseMatricesCSR

const DEFAULT_PARAMETER_FILE = joinpath(@__DIR__, "poisson_thyra.xml")

function main(distribute, parts)
    if !isfile(DEFAULT_PARAMETER_FILE)
        error(
            "Trilinos parameter XML file not found at $(DEFAULT_PARAMETER_FILE). " *
            "Add it to test/poisson_thyra.xml.",
        )
    end

    timer1 = MPI.Wtime()

    ranks = distribute(LinearIndices((prod(parts),)))
    domain = (0, 5, 0, 5)
    partition = (100, 20 * parts[2])
    model = CartesianDiscreteModel(ranks, parts, domain, partition)
    labels = get_face_labeling(model)
    add_tag_from_tags!(labels, "sides", [1, 2, 3, 4, 5, 6, 7, 8])

    order = 1
    reffe = ReferenceFE(lagrangian, Float64, order)
    V0 = TestFESpace(model, reffe; conformity=:H1, dirichlet_tags=["sides"])

    g(x) = 0.0
    Ug = TrialFESpace(V0, g)

    degree = 2
    Ω = Triangulation(model)
    dΩ = Measure(Ω, degree)

    f(x) = -1
    a(u, v) = ∫(∇(v) ⋅ ∇(u))dΩ
    l(v) = ∫(v * f)dΩ

    assem = SparseMatrixAssembler(SparseMatrixCSR{0,Float64,Int}, Vector{Float64}, Ug, V0)
    op = AffineFEOperator(a, l, Ug, V0, assem)

    timer2 = MPI.Wtime()
    solver = TrilinosSolve(DEFAULT_PARAMETER_FILE)
    uh = solve(solver, op)
    timer3 = MPI.Wtime()
    uh_lu = solve(LUSolver(), op)
    timer4 = MPI.Wtime()

    solution_difference_residual = sqrt(sum(∫(abs2(uh - uh_lu))dΩ))

    comm = MPI.COMM_WORLD
    rank = MPI.Comm_rank(comm)
    solving_time = MPI.Reduce(timer3 - timer2, MPI.MAX, 0, comm)
    lu_solving_time = MPI.Reduce(timer4 - timer3, MPI.MAX, 0, comm)
    gridap_setup_time = MPI.Reduce(timer2 - timer1, MPI.MAX, 0, comm)
    total_time = MPI.Reduce(timer4 - timer1, MPI.MAX, 0, comm)
    residual_value = solver.log.residual
    max_residual = MPI.Reduce(residual_value, MPI.MAX, 0, comm)
    max_solution_difference_residual = MPI.Reduce(solution_difference_residual, MPI.MAX, 0, comm)

    if rank == 0
        println("Trilinos solver time: ", solving_time)
        println("LU solver time: ", lu_solving_time)
        println("Gridap setup time: ", gridap_setup_time)
        println("Total program time: ", total_time)
        println("Solver: ", solver.log.name)
        println("Iterations: ", solver.log.num_iters)
        println("Residual: ", max_residual)
        println("Solution difference residual: ", max_solution_difference_residual)
    end

    return true
end

if abspath(PROGRAM_FILE) == @__FILE__
    with_mpi() do distribute
        rank_partition = (1, MPI.Comm_size(MPI.COMM_WORLD))
        main(distribute, rank_partition)
    end
end

# Run with MPI ranks:
# mpiexecjl --project=. -n 4 julia test/poisson_thyra.jl
