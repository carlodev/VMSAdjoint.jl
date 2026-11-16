# Tests for the chunked frozen-operator θ-method driver.
#
# The driver's job is bookkeeping: produce exactly the same sequence of steps,
# times and states as a plain `Gridap.solve(ThetaMethod, ...)`, while rebuilding
# the operator only every `freq` steps. On a problem whose forms really ARE
# constant, freezing is not an approximation, so agreement must be to machine
# precision for every chunk size — which is exactly what pins the step counting,
# the chunk-boundary state handoff and the time grid.

using Gridap
using Gridap.FESpaces
using VMSAdjoint.IncompressibleSolvers: march_frozen_theta

@testset "march_frozen_theta" begin
    model = CartesianDiscreteModel((0, 1, 0, 1), (4, 4))
    reffe = ReferenceFE(lagrangian, Float64, 1)
    V = TestFESpace(model, reffe; conformity=:H1, dirichlet_tags="boundary")
    U = TransientTrialFESpace(V, t -> (x -> 0.0))

    Ω = Triangulation(model)
    dΩ = Measure(Ω, 2)

    # heat equation: ∂u/∂t - Δu = f, constant in time
    m(t, u, v) = ∫(u * v)dΩ
    res(t, u, v) = ∫(∇(u) ⋅ ∇(v))dΩ
    rhs(t, v) = ∫(1.0 * v)dΩ

    t0, tF, dt, θ = 0.0, 0.6, 0.05, 0.5
    nsteps = round(Int, (tF - t0) / dt)

    u_init(x) = sin(pi * x[1]) * sin(pi * x[2])
    make_x0() = interpolate_everywhere(u_init, U(t0))

    # reference: one operator, no chunking, Gridap's own loop
    op_ref = TransientLinearFEOperator((res, m), rhs, U, V)
    sol_ref = Gridap.solve(ThetaMethod(LUSolver(), dt, θ), op_ref, t0, tF, make_x0())
    times_ref = Float64[]
    states_ref = Vector{Float64}[]
    for (t, uh) in sol_ref
        push!(times_ref, t)
        push!(states_ref, copy(get_free_dof_values(uh)))
    end
    @test length(times_ref) == nsteps

    make_op() = TransientLinearFEOperator((res, m), rhs, U, V; constant_forms=(true, true))

    for freq in (1, 2, 4, 7, 100)
        times = Float64[]
        states = Vector{Float64}[]
        steps = Int[]
        nrefresh = Ref(0)

        xF = march_frozen_theta(make_op, LUSolver(), U, make_x0(),
            t0, tF, dt, θ, freq;
            on_step! = (step, t, uh) -> begin
                push!(steps, step)
                push!(times, t)
                push!(states, copy(get_free_dof_values(uh)))
            end,
            refresh! = (_step, _state) -> (nrefresh[] += 1),
            label="test")

        # every step is produced exactly once, in order, on the same time grid
        @test steps == collect(1:nsteps)
        @test times ≈ times_ref
        # ... and the states match the unchunked reference to machine precision
        @test all(states[k] ≈ states_ref[k] for k in 1:nsteps)
        # the returned final state is the last one marched
        @test get_free_dof_values(xF) ≈ states_ref[end]
        # one operator rebuild per chunk
        @test nrefresh[] == cld(nsteps, freq)
    end

    # freq <= 0 is clamped to 1 rather than looping forever
    steps = Int[]
    march_frozen_theta(make_op, LUSolver(), U, make_x0(), t0, tF, dt, θ, 0;
        on_step! = (step, _t, _uh) -> push!(steps, step),
        refresh! = (_s, _v) -> nothing, label="test")
    @test steps == collect(1:nsteps)
end
