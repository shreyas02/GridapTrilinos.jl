using GridapTrilinos
using MPI
using PartitionedArrays
using Test

@testset "GridapTrilinos package" begin
    @test isdefined(GridapTrilinos, :TrilinosSolve)
    @test isdefined(GridapTrilinos, :ConstructTpetraMatrixWrapper)
    @test isdefined(GridapTrilinos, :ConstructTpetraVectorWrapper)
    @test isdefined(GridapTrilinos, :TrilinosSolveWrapper)
    @test isdefined(GridapTrilinos, :CopySolutionWrapper)
end

include("poisson_thyra.jl")
include("poisson_frosch.jl")
include("transient_cached.jl")

env_true(name) = lowercase(get(ENV, name, "false")) in ("1", "true", "yes")

@testset "Workflow setup" begin
    @test isdefined(Main, :main)
    @test DEFAULT_PARAMETER_FILE == joinpath(@__DIR__, "poisson_thyra.xml")
    @test isfile(DEFAULT_PARAMETER_FILE)
    @test isdefined(Main, :poisson_frosch)
    @test FROSCH_PARAMETER_FILE == joinpath(@__DIR__, "poisson_frosch.xml")
    @test isfile(FROSCH_PARAMETER_FILE)
    @test isdefined(Main, :transient_cached)
    @test TRANSIENT_PARAMETER_FILE == joinpath(@__DIR__, "poisson_frosch.xml")
    @test isfile(TRANSIENT_PARAMETER_FILE)
end

@testset "MPI workflows" begin
    workflows = (
        (name="Thyra", solve=main),
        (name="FROSch", solve=poisson_frosch),
        (name="Transient cached FROSch", solve=transient_cached),
    )

    if !env_true("GRIDAPTRILINOS_RUN_MPI_TESTS")
        @test_skip "Set GRIDAPTRILINOS_RUN_MPI_TESTS=true to run all MPI workflows."
    else
        with_mpi() do distribute
            comm_size = MPI.Comm_size(MPI.COMM_WORLD)
            rank_partition = (1, comm_size)

            for workflow in workflows
                @testset "$(workflow.name)" begin
                    result = workflow.solve(distribute, rank_partition)
                    GC.gc(true)
                    @test result
                end
            end
        end
    end
end
