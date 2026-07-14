using GridapTrilinos
using MPI
using PartitionedArrays
using Test

@testset "GridapTrilinos package" begin
    @test isdefined(GridapTrilinos, :TrilinosSolve)
    @test isdefined(GridapTrilinos, :TrilinosParallel)
end

include("poisson_frosch.jl")

@testset "Poisson FROSch workflow" begin
    @test isdefined(Main, :main)
    @test DEFAULT_PARAMETER_FILE == joinpath(@__DIR__, "poisson_frosch.xml")

    run_solver_test = lowercase(get(ENV, "GRIDAPTRILINOS_RUN_SOLVER_TEST", "false")) in ("1", "true", "yes")
    parameter_file = get(ENV, "GRIDAPTRILINOS_PARAMETER_FILE", DEFAULT_PARAMETER_FILE)
    residual_tol = parse(Float64, get(ENV, "GRIDAPTRILINOS_RESIDUAL_TOL", "1.0e-8"))

    if run_solver_test
        @test isfile(parameter_file)
        if isfile(parameter_file)
            with_mpi() do distribute
                parts = (2, 2)
                comm_size = MPI.Comm_size(MPI.COMM_WORLD)
                @test comm_size == prod(parts)
                if comm_size == prod(parts)
                    @test main(distribute, parts; parameter_file, residual_tol)
                end
            end
        end
    else
        @test_skip "Set GRIDAPTRILINOS_RUN_SOLVER_TEST=true and provide test/poisson_frosch.xml to run the MPI/Trilinos solve."
    end
end
