#############################################################################
# Adjoint incompressible Navier–Stokes solvers (steady / unsteady VMS)
#############################################################################

"""
    solve_inc_adj(am, simcase, d_bc, [filename,] timed::Symbol, uh, ph)

Solve the adjoint problem, `timed ∈ (:steady, :unsteady)`, with airfoil
Dirichlet value `d_bc = -dJ/d[CD,CL]` and primal solution `(uh, ph)`.
When `filename` is omitted it defaults to `"inc-adj-<timed>"`.
Returns the adjoint fields `(ϕu, ϕp)`; for `:unsteady` these are time-averaged,
and the time-averaged wall-shear correlation is stored in `am.params[:wall_shear_corr]`.
"""
function solve_inc_adj(am::AirfoilModel, simcase::Airfoil, d_bc::Vector{Float64}, timed::Symbol, uh, ph)
    solve_inc_adj(am, simcase, d_bc, "inc-adj-$timed", timed, uh, ph)
end

function solve_inc_adj(am::AirfoilModel, simcase::Airfoil, d_bc::Vector{Float64}, filename::String, timed::Symbol, uh, ph)
    if timed === :steady
        solve_inc_adj_steady(am, simcase, d_bc, filename, uh, ph)
    elseif timed === :unsteady
        solve_inc_adj_unsteady(am, simcase, d_bc, filename, uh, ph)
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
        dirichlet_masks=[(true, true), (true, true), (false, true)])
    reffe_p_adj = ReferenceFE(lagrangian, Float64, order)
    Q_adj = TestFESpace(model, reffe_p_adj, conformity=:H1, dirichlet_tags="inlet")

    return V_adj, Q_adj
end

"""
    register_boundary_measures!(params, model, order)

Register in `params` the measures and normals of the outer/airfoil boundaries
(`:dΓout/:nΓout`, `:dΓlim/:nΓlim`, `:dΓairfoil/:nΓairfoil`).
"""
function register_boundary_measures!(params::Dict{Symbol,Any}, model, order::Int)
    for (tag, suffix) in (("outlet", "out"), ("limits", "lim"), ("airfoil", "airfoil"))
        Γ = BoundaryTriangulation(model; tags=tag)
        updatekey(params, Symbol(:dΓ, suffix), Measure(Γ, 2 * order))
        updatekey(params, Symbol(:nΓ, suffix), get_normal_vector(Γ))
    end
end

#############################################################################
# Steady solver
#############################################################################

function solve_inc_adj_steady(am::AirfoilModel, simcase::Airfoil, d_boundary::Vector{Float64}, filename, uh, ph)
    @unpack params, model = am
    @unpack Ω = params
    @sunpack order = simcase

    register_boundary_measures!(params, model, order)

    V_adj, Q_adj = create_adjoint_spaces(model, simcase)
    @info "Solving Steady Adjoint, airfoil boundary condition: $d_boundary"

    U_adj = TrialFESpace(V_adj, [VectorValue(d_boundary...), VectorValue(0, 0), VectorValue(0, 0)])
    P_adj = TrialFESpace(Q_adj, 0.0)

    Y_adj = MultiFieldFESpace([V_adj, Q_adj])
    X_adj = MultiFieldFESpace([U_adj, P_adj])

    updatekey(params, :uh, uh)
    updatekey(params, :ph, ph)

    res_adj, rhs_adj = equations_adjoint(simcase, params, :steady)
    op_adj = AffineFEOperator(res_adj, rhs_adj, X_adj, Y_adj)

    solver = LinearFESolver(LUSolver())
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

function solve_inc_adj_unsteady(am::AirfoilModel, simcase::Airfoil, d_boundary::Vector{Float64}, filename, uh, ph)
    @unpack params, model = am
    @sunpack D, order, t0, tF, θ, dt, time_window = simcase
    @unpack Ω, UH = params

    register_boundary_measures!(params, model, order)

    V_adj, Q_adj = create_adjoint_spaces(model, simcase)
    @info "Solving Unsteady Adjoint, airfoil boundary condition: $d_boundary"

    u0 = ConstantInTime(VectorValue(d_boundary...))
    u_walls = ConstantInTime(VectorValue(zeros(D)...))
    p0 = ConstantInTime(0.0)

    U_adj = TransientTrialFESpace(V_adj, [u0, u_walls, u_walls])
    P_adj = TransientTrialFESpace(Q_adj, p0)

    Y_adj = MultiFieldFESpace([V_adj, Q_adj])
    X_adj = MultiFieldFESpace([U_adj, P_adj])

    uh0_adj = interpolate(u0(0.0), U_adj(0.0))
    ph0_adj = interpolate(p0(0.0), P_adj(0.0))

    # the adjoint marches backward: start from the LAST primal state
    copyto!(params[:uh].free_values, UH[end])

    xh0_adj = interpolate([uh0_adj, ph0_adj], X_adj(0.0))

    m_adj, res_adj, rhs_adj = equations_adjoint(simcase, params, :unsteady)
    op_adj = TransientLinearFEOperator((res_adj, m_adj), rhs_adj, X_adj, Y_adj)

    ode_solver = ThetaMethod(LUSolver(), dt, θ)
    sol = Gridap.solve(ode_solver, op_adj, t0, tF, xh0_adj)

    UH_ADJ = [copy(uh0_adj.free_values)]
    PH_ADJ = [copy(ph0_adj.free_values)]

    res_path = "Results_unsteady_primal"
    mkpath(res_path)

    time_vec = collect(t0:dt:tF)
    t_length = length(time_vec)

    createpvd(filename) do pvd
        pvd[time_vec[end]] = createvtk(Ω, joinpath(res_path, "$(filename)_$tF.vtu");
            nsubcells=order, cellfields=["uh-adj" => uh0_adj, "ph-adj" => ph0_adj])
        for (idx, (t, xhtn)) in enumerate(sol)
            ϕu, ϕp = xhtn

            push!(UH_ADJ, copy(ϕu.free_values))
            push!(PH_ADJ, copy(ϕp.free_values))

            # adjoint pseudo-time t maps to physical time t_adj = tF - t
            idx_adj = t_length - idx
            t_adj = time_vec[idx_adj]

            # feed the primal state at the SAME physical time into the adjoint operator
            idx_adj > 0 && copyto!(params[:uh].free_values, UH[idx_adj])

            println("Adjoint solved at time step $t_adj")

            pvd[t_adj] = createvtk(Ω, joinpath(res_path, "$(filename)_$t_adj.vtu");
                nsubcells=order, cellfields=["uh-adj" => ϕu, "ph-adj" => ϕp])
        end
    end

    jldsave("UnsteadyAdjointFields.jld2"; UH_ADJ, PH_ADJ)

    # average over the (reversed) physical window, skipping the adjoint start-up
    time_window_adj = (tF - time_window[1], tF - 10 * dt)
    @assert time_window_adj[2] > time_window_adj[1] "Adjoint Time Window Averaging not consistent"

    avg_UH_ADJ, avg_PH_ADJ = time_average_fields(UH_ADJ, PH_ADJ, time_window_adj, dt, t0)
    uh0_adj.free_values .= avg_UH_ADJ
    ph0_adj.free_values .= avg_PH_ADJ

    store_wall_shear_correlation!(am, simcase, U_adj(0.0), UH, UH_ADJ, time_vec)

    return uh0_adj, ph0_adj
end

"""
    store_wall_shear_correlation!(am, simcase, Uadj0, UH, UH_ADJ, time_vec)

Covariance-correct unsteady sensitivity seed (Srinath & Mittal, JCP 2010).

Builds the time-AVERAGED wall-shear product `S = ⟨(∂v_t/∂n)(∂ψ_t/∂n)⟩` over the
physical averaging window, pairing primal `v(t)` and adjoint `ψ(t)` at the SAME
physical time, and stores it in `am.params[:wall_shear_corr]`. Storing `S`
(rather than multiplying the two time-averaged fields) keeps the primal–adjoint
covariance term; `compute_gradient` then only needs `∫_Γ ν S δβ dΓ` per design
parameter.

The primal is rebuilt on the test space (airfoil no-slip ⇒ homogeneous Dirichlet),
the adjoint on the trial space at t=0 (airfoil Dirichlet constant in time).
"""
function store_wall_shear_correlation!(am::AirfoilModel, simcase::Airfoil, Uadj0, UH, UH_ADJ, time_vec)
    @unpack params, model = am
    @unpack nΓ = params
    @sunpack tF, dt, time_window = simcase

    V_prim, _ = create_primal_spaces(model, simcase)
    ∂ₙᵗ(u) = tangential_normal_derivative(u, nΓ)

    # skip the last 10*dt near tF: that is the adjoint start-up transient
    t_hi = min(time_window[2], tF - 10 * dt)
    win = findall(t -> time_window[1] <= t <= t_hi, time_vec)
    @assert !isempty(win) "Sensitivity averaging window contains no time step"

    # physical index j: primal = UH[j], adjoint = UH_ADJ[t_length - j + 1]
    t_length = length(time_vec)
    wall_shear_products = map(win) do j
        v_j = FEFunction(V_prim, UH[j])
        ψ_j = FEFunction(Uadj0, UH_ADJ[t_length-j+1])
        ∂ₙᵗ(v_j) ⋅ ∂ₙᵗ(ψ_j)
    end

    updatekey(params, :wall_shear_corr, sum(wall_shear_products) / length(wall_shear_products))
end
