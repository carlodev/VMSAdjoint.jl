using VMSAdjoint
using AirfoilTools
using Gridap
using Gridap.Algebra
using SparseArrays
using Test

using GridapSolvers
using VMSAdjoint.IncompressibleSolvers: ConstantInTime, time_average_fields,
    full_mask, wall_normal_mask, boundary_vector, set_free_values!,
    linear_solver, GPULUSolver, SolverContext, schur_biform
using VMSAdjoint.IteratorTools: bounds_w

@testset "VMSAdjoint.jl" begin

    @testset "ConstantInTime" begin
        f = ConstantInTime(VectorValue(1.0, 2.0))
        @test f(VectorValue(0.3, 0.4), 12.5) == VectorValue(1.0, 2.0)
        @test f(12.5)(VectorValue(0.3, 0.4)) == VectorValue(1.0, 2.0)
    end

    @testset "Dirichlet helpers" begin
        @test full_mask(2) == (true, true)
        @test full_mask(3) == (true, true, true)
        @test wall_normal_mask(2) == (false, true)
        @test wall_normal_mask(3) == (false, true, false)
        # 2-component adjoint BC zero-padded to the spanwise direction
        @test boundary_vector([1.5, -2.0], 2) == VectorValue(1.5, -2.0)
        @test boundary_vector([1.5, -2.0], 3) == VectorValue(1.5, -2.0, 0.0)
    end

    @testset "time_average_fields" begin
        # 5 snapshots at t = 0, 0.1, ..., 0.4
        UH = [fill(Float64(k), 3) for k in 0:4]
        PH = [fill(Float64(-k), 2) for k in 0:4]
        avg_u, avg_p = time_average_fields(UH, PH, (0.1, 0.3), 0.1, 0.0)
        @test avg_u ≈ fill(2.0, 3)  # mean of snapshots 1,2,3
        @test avg_p ≈ fill(-2.0, 2)
        # window ↔ snapshot mapping honours t0
        avg_u, _ = time_average_fields(UH, PH, (1.1, 1.2), 0.1, 1.0)
        @test avg_u ≈ fill(1.5, 3)
    end

    @testset "linear solver backends" begin
        @test linear_solver(:cpu) isa LUSolver
        @test linear_solver(:gpu) isa GPULUSolver
        @test_throws ArgumentError linear_solver(:tpu)
        # without CUDA/CUDSS loaded the GPU backend must fail loudly, not silently
        A = sparse([1.0 0.0; 0.0 1.0])
        @test_throws ErrorException Gridap.Algebra.symbolic_setup(GPULUSolver(), A)
        @test AdjSolver(linear_backend=:gpu).linear_backend === :gpu
        @test AdjSolver(linear_backend=:petsc).linear_backend === :petsc
        @test AdjSolver(linear_backend=:block).linear_backend === :block
        @test_throws Exception AdjSolver(linear_backend=:apu)

        # :block needs the pressure spaces; asking for it without them must say so
        @test_throws ErrorException linear_solver(:block)

        model = CartesianDiscreteModel((0, 1, 0, 1), (2, 2))
        Q = TestFESpace(model, ReferenceFE(lagrangian, Float64, 1);
            conformity=:H1, dirichlet_tags="boundary")
        P = TrialFESpace(Q, 0.0)
        dΩ = Measure(Triangulation(model), 2)
        ctx = SolverContext(P=P, Q=Q, dΩ=dΩ, ν=0.01, τ=0.05)
        @test linear_solver(:block, ctx) isa GridapSolvers.LinearSolvers.FGMRESSolver

        # the Schur approximation interpolates between the two limits:
        # τ = θ·dt weights the pressure Laplacian, 1/ν the pressure mass
        s_transient = schur_biform(ctx)
        s_steady = schur_biform(SolverContext(P=P, Q=Q, dΩ=dΩ, ν=0.01, τ=Inf))
        p = interpolate(x -> x[1], P)
        # steady drops the Laplacian term entirely (no mass matrix in the operator)
        @test sum(s_steady(p, p)) < sum(s_transient(p, p))
        @test isfinite(sum(s_steady(p, p)))
    end

    @testset "mesh motion options" begin
        @test AdjSolver().mesh_motion === true
        @test AdjSolver(mesh_motion=false).mesh_motion === false
        @test AdjSolver(mesh_stiffening=2.0).mesh_stiffening == 2.0
        @test_throws Exception AdjSolver(mesh_stiffening=-1.0)
    end

    @testset "CST bounds honour DesignBounds" begin
        Nhalf = 5
        cstd = AirfoilCSTDesign(CSTGeometry(cstw=CSTweights(Nhalf, 0.2)), 40)
        bounds = DesignBounds(upper=1.5, lower=-1.5, nose=0.08, cross=0.15)
        lb, ub = bounds_w(cstd, 2 * Nhalf, bounds)
        @test length(lb) == length(ub) == 2 * Nhalf
        @test all(ub .> lb)
        @test lb[1:3] == fill(0.08, 3)             # leading-edge suction weights >= nose
        @test lb[4:Nhalf] == fill(-0.15, Nhalf - 3) # remaining suction weights >= -cross
        @test all(lb[Nhalf+1:end] .== -1.5)         # pressure side >= lower
        @test all(ub[1:Nhalf] .== 1.5)              # suction side <= upper
        # user values are actually consumed
        lb2, ub2 = bounds_w(cstd, 2 * Nhalf, DesignBounds(upper=0.7, lower=-0.7))
        @test all(ub2[1:Nhalf] .== 0.7) && all(lb2[Nhalf+1:end] .== -0.7)
    end

    @testset "mesh motion" begin
        include("mesh_motion.jl")
    end

    @testset "transient solve" begin
        include("transient_solve.jl")
    end

    # Heavy integration smoke test (mesh generation + steady primal/adjoint/sensitivity).
    # Run with: VMSADJOINT_SMOKE=1 julia --project -e 'using Pkg; Pkg.test()'
    if get(ENV, "VMSADJOINT_SMOKE", "0") == "1"
        @testset "steady smoke" begin
            include("smoke_steady.jl")
        end
    end

end
