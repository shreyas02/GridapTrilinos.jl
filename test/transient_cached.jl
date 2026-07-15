using Gridap
using GridapDistributed
using GridapTrilinos
using MPI
using PartitionedArrays
using SparseMatricesCSR

const TRANSIENT_PARAMETER_FILE = joinpath(@__DIR__, "poisson_frosch.xml")

function transient_cached(distribute, parts; t0=0.0, tF=0.20, Δt=0.05, θ=0.5)
    if !isfile(TRANSIENT_PARAMETER_FILE)
        error(
            "Trilinos parameter XML file not found at $(TRANSIENT_PARAMETER_FILE). " *
            "Add it to test/poisson_frosch.xml.",
        )
    end

    ranks = distribute(LinearIndices((prod(parts),)))
    domain = (-1, +1, -1, +1)
    partition = (20, 20 * parts[2])
    model = CartesianDiscreteModel(ranks, parts, domain, partition)

    order = 1
    reffe = ReferenceFE(lagrangian, Float64, order)
    V0 = TestFESpace(model, reffe; conformity=:H1, dirichlet_tags="boundary")

    g(t) = x -> exp(-2 * t) * sinpi(t * x[1]) * (x[2]^2 - 1)
    Ug = TransientTrialFESpace(V0, g)

    degree = 2
    Ω = Triangulation(model)
    dΩ = Measure(Ω, degree)

    α(t) = x -> 1.0
    f(t) = x -> sin(t) * sinpi(x[1]) * sinpi(x[2])

    m(t, dtu, v) = ∫(v * dtu)dΩ
    a(t, u, v) = ∫(α(t) * ∇(v) ⋅ ∇(u))dΩ
    l(t, v) = ∫(v * f(t))dΩ

    assem = SparseMatrixAssembler(SparseMatrixCSR{0,Float64,Int}, Vector{Float64}, Ug, V0)
    op = TransientLinearFEOperator(
        (a, m),
        l,
        Ug,
        V0,
        constant_forms=(true, true),
        assembler=assem,
    )

    trilinos_linear_solver = TrilinosSolve(TRANSIENT_PARAMETER_FILE)
    trilinos_solver = ThetaMethod(trilinos_linear_solver, Δt, θ)

    uh0 = interpolate_everywhere(g(t0), Ug(t0))
    uh_trilinos = solve(trilinos_solver, op, t0, tF, uh0)
    uh_lu = solve(ThetaMethod(LUSolver(), Δt, θ), op, t0, tF, uh0)

    final_time = t0
    trilinos_final_energy = 0.0
    lu_final_energy = 0.0
    final_difference = 0.0
    num_steps = 0
    for ((tn, uhn), (tn_lu, uhn_lu)) in zip(uh_trilinos, uh_lu)
        @assert isapprox(tn, tn_lu)
        final_time = tn
        trilinos_final_energy = sqrt(sum(∫(abs2(uhn))dΩ))
        lu_final_energy = sqrt(sum(∫(abs2(uhn_lu))dΩ))
        final_difference = sqrt(sum(∫(abs2(uhn - uhn_lu))dΩ))
        num_steps += 1
    end

    comm = MPI.COMM_WORLD
    rank = MPI.Comm_rank(comm)
    max_trilinos_final_energy = MPI.Reduce(trilinos_final_energy, MPI.MAX, 0, comm)
    max_lu_final_energy = MPI.Reduce(lu_final_energy, MPI.MAX, 0, comm)
    max_final_difference = MPI.Reduce(final_difference, MPI.MAX, 0, comm)
    max_residual = MPI.Reduce(trilinos_linear_solver.log.residual, MPI.MAX, 0, comm)

    if rank == 0
        println("Transient final time: ", final_time)
        println("Transient steps: ", num_steps)
        println("Solver: ", trilinos_linear_solver.log.name)
        println("Last solve iterations: ", trilinos_linear_solver.log.num_iters)
        println("Last solve residual: ", max_residual)
        println("Trilinos final solution L2 norm: ", max_trilinos_final_energy)
        println("LU final solution L2 norm: ", max_lu_final_energy)
        println("Final Trilinos-LU L2 difference: ", max_final_difference)
    end

    uh_trilinos = nothing
    uh_lu = nothing
    uh0 = nothing
    trilinos_solver = nothing
    trilinos_linear_solver = nothing
    op = nothing
    GC.gc(true)

    return num_steps > 0 && (rank != 0 || max_final_difference < 1.0e-6)
end

if abspath(PROGRAM_FILE) == @__FILE__
    with_mpi() do distribute
        rank_partition = (1, MPI.Comm_size(MPI.COMM_WORLD))
        transient_cached(distribute, rank_partition)
        GC.gc(true)
    end
end

# Run with MPI ranks:
# mpiexecjl --project=. -n 4 julia test/transient_cached.jl
