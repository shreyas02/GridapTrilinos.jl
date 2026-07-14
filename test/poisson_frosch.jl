using Gridap
using GridapDistributed
using GridapTrilinos
using MPI
using PartitionedArrays

const DEFAULT_PARAMETER_FILE = joinpath(@__DIR__, "poisson_frosch.xml")

function main(distribute, parts; parameter_file=DEFAULT_PARAMETER_FILE, residual_tol=1.0e-8)
    if !isfile(parameter_file)
        error(
            "Trilinos parameter XML file not found at $(parameter_file). " *
            "Pass the XML path as the first command-line argument or add it to test/poisson_frosch.xml.",
        )
    end

    timer1 = MPI.Wtime()

    ranks = distribute(LinearIndices((prod(parts),)))
    domain = (0, 5, 0, 5)
    partition = (100, 100)
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

    op = AffineFEOperator(a, l, Ug, V0)

    timer2 = MPI.Wtime()
    solver = TrilinosSolve(parameter_file)
    uh = solve(solver, op)
    timer3 = MPI.Wtime()

    comm = MPI.COMM_WORLD
    rank = MPI.Comm_rank(comm)
    solving_time = MPI.Reduce(timer3 - timer2, MPI.MAX, 0, comm)
    gridap_setup_time = MPI.Reduce(timer2 - timer1, MPI.MAX, 0, comm)
    total_time = MPI.Reduce(timer3 - timer1, MPI.MAX, 0, comm)
    residual_value = solver.log.residual
    max_residual = MPI.Reduce(residual_value, MPI.MAX, 0, comm)

    if rank == 0
        println("Trilinos solver time: ", solving_time)
        println("Gridap setup time: ", gridap_setup_time)
        println("Total program time: ", total_time)
        println("Solver: ", solver.log.name)
        println("Iterations: ", solver.log.num_iters)
        println("Residual: ", max_residual)
    end

    if rank == 0
        return isfinite(max_residual) && max_residual <= residual_tol
    end
    return true
end

if abspath(PROGRAM_FILE) == @__FILE__
    rank_partition = (2, 2)
    parameter_file = isempty(ARGS) ? DEFAULT_PARAMETER_FILE : ARGS[1]

    with_mpi() do distribute
        main(distribute, rank_partition; parameter_file)
    end
end

# Run with 4 MPI ranks:
# mpiexecjl --project=. -n 4 julia test/poisson_frosch.jl
#
# To use a different XML file:
# mpiexecjl --project=. -n 4 julia test/poisson_frosch.jl test/my_parameters.xml
