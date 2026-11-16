#############################################################################
# Adjoint incompressible Navier–Stokes solvers (steady / unsteady VMS)
#############################################################################

"""
    solve_inc_adj(am, simcase, d_bc, [filename,] timed::Symbol, uh, ph;
                  linear_backend=:cpu, adj_transient_steps=10)

Solve the adjoint problem, `timed ∈ (:steady, :unsteady)`, with airfoil
Dirichlet value `d_bc = -dJ/d[CD,CL]` and primal solution `(uh, ph)`.
When `filename` is omitted it defaults to `"inc-adj-<timed>"`.
`linear_backend ∈ (:cpu, :gpu)` selects the linear solver (see [`linear_solver`](@ref)).
Returns the adjoint fields `(ϕu, ϕp)`; for `:unsteady` these are time-averaged,
and the time-averaged wall-shear correlation is stored in `am.params[:wall_shear_corr]`.

`adj_transient_steps` is how many time steps of adjoint start-up transient are
discarded from the averaging (see [`solve_inc_adj_unsteady`](@ref)); it is
ignored by the steady solver.
"""
function solve_inc_adj(am::AirfoilModel, simcase::Airfoil, d_bc::Vector{Float64}, timed::Symbol, uh, ph; linear_backend::Symbol=:cpu, adj_transient_steps::Int=10)
    solve_inc_adj(am, simcase, d_bc, "inc-adj-$timed", timed, uh, ph; linear_backend=linear_backend, adj_transient_steps=adj_transient_steps)
end

function solve_inc_adj(am::AirfoilModel, simcase::Airfoil, d_bc::Vector{Float64}, filename::String, timed::Symbol, uh, ph; linear_backend::Symbol=:cpu, adj_transient_steps::Int=10)
    if timed === :steady
        solve_inc_adj_steady(am, simcase, d_bc, filename, uh, ph, linear_backend)
    elseif timed === :unsteady
        solve_inc_adj_unsteady(am, simcase, d_bc, filename, uh, ph, linear_backend, adj_transient_steps)
    else
        throw(ArgumentError("timed must be :steady or :unsteady, got :$timed"))
    end
end

#############################################################################
# Spaces and measures
#############################################################################

function create_adjoint_spaces(model, simcase::Airfoil)
    @sunpack order, D = simcase
    reffe_u_adj = ReferenceFE(lagrangian, VectorValue{D,Float64}, order)
    V_adj = TestFESpace(model, reffe_u_adj, conformity=:H1,
        dirichlet_tags=["airfoil", "outlet", "limits"],
        dirichlet_masks=[full_mask(D), full_mask(D), wall_normal_mask(D)])
    reffe_p_adj = ReferenceFE(lagrangian, Float64, order)
    Q_adj = TestFESpace(model, reffe_p_adj, conformity=:H1, dirichlet_tags="inlet")

    return V_adj, Q_adj
end

#############################################################################
# Steady solver
#############################################################################

function solve_inc_adj_steady(am::AirfoilModel, simcase::Airfoil, d_boundary::Vector{Float64}, filename, uh, ph, linear_backend::Symbol=:cpu)
    @unpack params, model = am
    @unpack Ω = params
    @sunpack order, D = simcase

    V_adj, Q_adj = create_adjoint_spaces(model, simcase)
    @info "Solving Steady Adjoint, airfoil boundary condition: $d_boundary"

    u_zero = VectorValue(zeros(D)...)
    U_adj = TrialFESpace(V_adj, [boundary_vector(d_boundary, D), u_zero, u_zero])
    P_adj = TrialFESpace(Q_adj, 0.0)

    Y_adj = velocity_pressure_space(V_adj, Q_adj, linear_backend)
    X_adj = velocity_pressure_space(U_adj, P_adj, linear_backend)

    updatekey(params, :uh, uh)
    updatekey(params, :ph, ph)

    res_adj, rhs_adj = equations_adjoint(simcase, params, :steady)
    op_adj = AffineFEOperator(res_adj, rhs_adj, X_adj, Y_adj)

    ctx = solver_context(simcase, params, P_adj, Q_adj; transient=false)
    solver = LinearFESolver(linear_solver(linear_backend, ctx))
    ϕu, ϕp = Gridap.solve(solver, op_adj)

    res_path = "Results_adj"
    mkpath(res_path)

    if !isnothing(filename)
        writevtk(Ω, joinpath(res_path, "$filename.vtu"); nsubcells=order,
            cellfields=["phi-u" => ϕu, "phi-p" => ϕp, "uh" => uh, "ph" => ph])
    end

    return ϕu, ϕp
end

#############################################################################
# Unsteady solver
#############################################################################

"""
    solve_inc_adj_unsteady(am, simcase, d_boundary, filename, uh, ph,
                           linear_backend=:cpu, adj_transient_steps=10)

Unsteady adjoint solve. The adjoint marches BACKWARD in physical time, so its
start-up transient sits at the *end* of the physical record, in
`t ∈ (tF - adj_transient_steps·dt, tF]`. Those steps are excluded both from the
time-average of the adjoint fields and from the wall-shear correlation, so that
the two are averaged over the same physical window.

`adj_transient_steps` is a convergence parameter, not a physical constant: how
long the adjoint takes to forget its terminal condition depends on the flow and
the discretization. Sweep it (together with the primal averaging window) before
trusting an unsteady gradient — a gradient that moves when it changes is
windowed rather than converged.
"""
function solve_inc_adj_unsteady(am::AirfoilModel, simcase::Airfoil, d_boundary::Vector{Float64}, filename, uh, ph, linear_backend::Symbol=:cpu, adj_transient_steps::Int=10)
    adj_transient_steps >= 0 || throw(ArgumentError(
        "adj_transient_steps must be non-negative, got $adj_transient_steps"))
    @unpack params, model = am
    @sunpack D, order, t0, tF, θ, dt, time_window = simcase
    @sunpack M = simcase # save a .vtu file every M steps
    @sunpack matrix_freq_update = simcase # rebuild+refactorize the operator every N steps
    @unpack Ω, UH = params

    V_adj, Q_adj = create_adjoint_spaces(model, simcase)
    @info "Solving Unsteady Adjoint, airfoil boundary condition: $d_boundary"

    u0 = ConstantInTime(boundary_vector(d_boundary, D))
    u_walls = ConstantInTime(VectorValue(zeros(D)...))
    p0 = ConstantInTime(0.0)

    U_adj = TransientTrialFESpace(V_adj, [u0, u_walls, u_walls])
    P_adj = TransientTrialFESpace(Q_adj, p0)

    Y_adj = velocity_pressure_space(V_adj, Q_adj, linear_backend)
    X_adj = velocity_pressure_space(U_adj, P_adj, linear_backend)

    uh0_adj = interpolate(u0(t0), U_adj(t0))
    ph0_adj = interpolate(p0(t0), P_adj(t0))

    xh0_adj = interpolate([uh0_adj, ph0_adj], X_adj(t0))

    UH_ADJ = [copy(uh0_adj.free_values)]
    PH_ADJ = [copy(ph0_adj.free_values)]

    res_path = "Results_unsteady_primal"
    mkpath(res_path)

    time_vec = collect(t0:dt:tF)
    t_length = length(time_vec)

    # Frozen linearization for the adjoint.
    #
    # Adjoint pseudo-step `idx` corresponds to physical index `t_length - idx`.
    # Unlike the primal — which can only lag its linearization, being causal —
    # the adjoint has the whole primal history available, so the frozen state is
    # taken at the MIDPOINT of the chunk it will be used for rather than at its
    # start. That centres the O(freq·dt) linearization error instead of biasing
    # it to one side, at no cost.
    #
    # NOTE: freezing the primal state over a chunk is a physics approximation
    # that feeds straight into the sensitivity. Validate a production
    # `matrix_freq_update` against `matrix_freq_update = 1` on one design before
    # trusting the gradient (see `march_frozen_theta`).
    set_linearization! = function (steps_done::Int)
        nnext = min(matrix_freq_update, (t_length - 1) - steps_done)
        nnext <= 0 && return nothing
        idx_mid = t_length - (steps_done + (nnext + 1) ÷ 2)
        idx_mid = clamp(idx_mid, 1, length(UH))
        copyto!(params[:uh].free_values, UH[idx_mid])
        return nothing
    end
    set_linearization!(0)

    ctx = solver_context(simcase, params, P_adj, Q_adj; transient=true)
    lsolver = linear_solver(linear_backend, ctx)

    make_operator = function ()
        m_adj, res_adj, rhs_adj = equations_adjoint(simcase, params, :unsteady)
        return TransientLinearFEOperator((res_adj, m_adj), rhs_adj, X_adj, Y_adj;
            constant_forms=(true, true))
    end

    refresh! = function (steps_done, _state)
        set_linearization!(steps_done)
    end

    createpvd(filename) do pvd
        pvd[time_vec[end]] = createvtk(Ω, joinpath(res_path, "$(filename)_$tF.vtu");
            nsubcells=order, cellfields=["uh-adj" => uh0_adj, "ph-adj" => ph0_adj])

        on_step! = function (idx, _t, xhtn)
            ϕu, ϕp = xhtn

            push!(UH_ADJ, copy(ϕu.free_values))
            push!(PH_ADJ, copy(ϕp.free_values))

            if mod(idx, M) == 0
                # adjoint pseudo-step idx maps to physical time time_vec[t_length - idx]
                t_adj = time_vec[max(1, t_length - idx)]
                pvd[t_adj] = createvtk(Ω, joinpath(res_path, "$(filename)_$t_adj.vtu");
                    nsubcells=order, cellfields=["uh-adj" => ϕu, "ph-adj" => ϕp])
            end
        end

        march_frozen_theta(make_operator, lsolver, X_adj, xh0_adj, t0, tF, dt, θ, matrix_freq_update;
            on_step! = on_step!, refresh! = refresh!, label = "Adjoint")
    end

    jldsave("UnsteadyAdjointFields.jld2"; UH_ADJ, PH_ADJ)

    # UH_ADJ is indexed by the adjoint pseudo-time τ ∈ [t0, tF]; physical time is
    # t = tF - (τ - t0). Average over the primal averaging window in physical time,
    # skipping the adjoint start-up transient (physical times in
    # (tF - adj_transient_steps·dt, tF], i.e. the first `adj_transient_steps`
    # pseudo-time steps) — same window as `store_wall_shear_correlation!`.
    t_hi_phys = min(time_window[2], tF - adj_transient_steps * dt)
    @assert t_hi_phys > time_window[1] "Adjoint Time Window Averaging not consistent: " *
                                       "time_window = $time_window, tF = $tF, adj_transient_steps = $adj_transient_steps"
    time_window_adj = (t0 + tF - t_hi_phys, t0 + tF - time_window[1])

    avg_UH_ADJ, avg_PH_ADJ = time_average_fields(UH_ADJ, PH_ADJ, time_window_adj, dt, t0)
    uh0_adj.free_values .= avg_UH_ADJ
    ph0_adj.free_values .= avg_PH_ADJ

    store_wall_shear_correlation!(am, simcase, U_adj(0.0), UH, UH_ADJ, time_vec, adj_transient_steps)

    return uh0_adj, ph0_adj
end

"""
    store_wall_shear_correlation!(am, simcase, Uadj0, UH, UH_ADJ, time_vec, adj_transient_steps=10)

Covariance-correct unsteady sensitivity seed (Srinath & Mittal, JCP 2010).

Builds the time-AVERAGED wall-shear product `S = ⟨(∂v_t/∂n)(∂ψ_t/∂n)⟩` over the
physical averaging window, pairing primal `v(t)` and adjoint `ψ(t)` at the SAME
physical time, and stores it in `am.params[:wall_shear_corr]`. Storing `S`
(rather than multiplying the two time-averaged fields) keeps the primal–adjoint
covariance term; `compute_gradient` then only needs `∫_Γ ν S δβ dΓ` per design
parameter.

The primal is rebuilt on the test space (airfoil no-slip ⇒ homogeneous Dirichlet),
the adjoint on the trial space at t=0 (airfoil Dirichlet constant in time).

`S` is accumulated as a Galerkin projection onto a discontinuous space on `Γ`,
NOT as `sum` over a vector of lazy `CellField`s. The latter builds an `N`-deep
nested `Operation` tree — every one of the `N-1` additions creates a deeper
parametric type that Julia specializes from scratch — and the tree is then
re-evaluated once per design parameter inside `sensitivity_loop`. At the window
lengths this paper uses (`N ≈ 500`) that does not terminate in practical time.
"""
function store_wall_shear_correlation!(am::AirfoilModel, simcase::Airfoil, Uadj0, UH, UH_ADJ, time_vec, adj_transient_steps::Int=10)
    @unpack params, model = am
    @unpack nΓ, Γ, dΓ = params
    @sunpack order, tF, dt, time_window = simcase

    V_prim, _ = create_primal_spaces(model, simcase)
    ∂ₙᵗ(u) = tangential_normal_derivative(u, nΓ)

    # skip the last adj_transient_steps*dt near tF: that is the adjoint start-up transient
    t_hi = min(time_window[2], tF - adj_transient_steps * dt)
    win = findall(t -> time_window[1] <= t <= t_hi, time_vec)
    @assert !isempty(win) "Sensitivity averaging window contains no time step"

    # physical index j: primal = UH[j], adjoint = UH_ADJ[t_length - j + 1]
    t_length = length(time_vec)
    S_snapshot(j) = ∂ₙᵗ(FEFunction(V_prim, UH[j])) ⋅
                    ∂ₙᵗ(FEFunction(Uadj0, UH_ADJ[t_length-j+1]))

    # On straight-sided faces ∂ₙᵗ of an order-`order` velocity is a polynomial of
    # degree `order-1`, so S has degree 2(order-1) per face and lies exactly in
    # `V_S`: the projection is an identity on S, not an approximation. `dΓ`
    # integrates with degree 2·order, so mass matrix and right-hand sides are
    # exact too. The result is a plain FEFunction — depth 1, one assembly per
    # snapshot, O(1) per design parameter afterwards.
    deg_S = max(2 * (order - 1), 0)
    V_S = FESpace(Γ, ReferenceFE(lagrangian, Float64, deg_S); conformity=:L2)
    M_fact = lu(assemble_matrix((u, v) -> ∫(u * v)dΓ, V_S, V_S))

    # Verify exactness on a few snapshots before spending the whole window on it:
    # if these are not at roundoff the degree argument is wrong for this mesh
    # (curved boundary faces) and every gradient built on S would be silently off.
    l2Γ(f) = sqrt(abs(sum(∫(f * f)dΓ)))
    for j in win[unique([1, cld(length(win), 2), length(win)])]
        Sx = S_snapshot(j)
        Sp = FEFunction(V_S, M_fact \ assemble_vector(v -> ∫(Sx * v)dΓ, V_S))
        err = l2Γ(Sp - Sx) / l2Γ(Sx)
        err < 1e-10 || error("wall-shear projection is not exact (relative L2 error " *
                             "$(err) at t = $(time_vec[j])): S is not in the " *
                             "degree-$(deg_S) space on this mesh.")
    end

    b_S = zeros(num_free_dofs(V_S))
    for j in win
        b_S .+= assemble_vector(v -> ∫(S_snapshot(j) * v)dΓ, V_S)
    end

    updatekey(params, :wall_shear_corr, FEFunction(V_S, M_fact \ (b_S ./ length(win))))
end
