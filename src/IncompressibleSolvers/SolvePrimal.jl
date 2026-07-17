#############################################################################
# Primal incompressible Navier–Stokes solvers (steady / unsteady VMS)
#############################################################################

"""
    ConstantInTime(value)

Callable representing a time-independent Dirichlet value, in the two-form
convention expected by Gridap transient spaces: `f(x, t)` and `f(t::Real) -> x -> value`.
Used for inlet/wall boundary conditions of both the primal and the adjoint problem.
"""
struct ConstantInTime{T} <: Function
    value::T
end
(f::ConstantInTime)(x, t) = f.value
(f::ConstantInTime)(t::Real) = x -> f.value

"""
    tangential_normal_derivative(u, nΓ)

Tangential projection of the wall-normal derivative of the vector field `u` on a
boundary with outward normal `nΓ`:

    ∂ₙᵗu = (I - n ⊗ n) ⋅ (∂u/∂n),    (∂u/∂n)ᵢ = nⱼ ∂uᵢ/∂xⱼ

In Gridap's convention (∇u)[i,j] = ∂u_j/∂x_i, the normal derivative is
`transpose(∇(u)) ⋅ nΓ`; subtracting its normal component leaves the tangential
shear vector.

Replaces the former 2D-only kernel

    ∂ₜ∂n(u) = (∇(u) ⋅ tΓ) ⋅ nΓ,      tΓ = nΓ rotated by +90°

which returned the scalar shear ∂u_t/∂n along the single in-plane tangent `tΓ`.
The two formulations produce identical sensitivity kernels in 2D:

    ∂ₜ∂n(u) * ∂ₜ∂n(ψ) == ∂ₙᵗu ⋅ ∂ₙᵗψ

since in 2D the tangential space is spanned by `tΓ` alone (and at a no-slip wall
the normal-normal component ∂u_n/∂n vanishes by incompressibility — the projector
merely removes its discretization residue). The projector form needs no tangent
vector, so it carries over unchanged to 3D, where a single rotated tangent would
drop the spanwise shear contribution ν(∂u_z/∂n)(∂ψ_z/∂n).
"""
function tangential_normal_derivative(u, nΓ)
    ∂ₙu = transpose(∇(u)) ⋅ nΓ
    return ∂ₙu - (∂ₙu ⋅ nΓ) * nΓ
end

#############################################################################
# Entry points
#############################################################################

"""
    solve_inc_primal(am, simcase, [filename,] timed::Symbol; uh0=nothing, ph0=nothing)

Solve the primal incompressible flow problem, `timed ∈ (:steady, :unsteady)`.
When `filename` is omitted it defaults to `"inc-primal-<timed>"`.
Returns `(uh, ph)`; for `:unsteady` these are the time-averaged fields over
`time_window`, and the full history is stored in `am.params[:UH]/[:PH]`.
"""
function solve_inc_primal(am::AirfoilModel, simcase::Airfoil, timed::Symbol; uh0=nothing, ph0=nothing)
    solve_inc_primal(am, simcase, "inc-primal-$timed", timed; uh0=uh0, ph0=ph0)
end

function solve_inc_primal(am::AirfoilModel, simcase::Airfoil, filename::String, timed::Symbol; uh0=nothing, ph0=nothing)
    if timed === :steady
        solve_inc_primal_steady(am, simcase, filename, uh0, ph0)
    elseif timed === :unsteady
        solve_inc_primal_unsteady(am, simcase, filename, uh0, ph0)
    else
        throw(ArgumentError("timed must be :steady or :unsteady, got :$timed"))
    end
end

#############################################################################
# Spaces and measures
#############################################################################

function create_primal_spaces(model, simcase::Airfoil)
    @sunpack order, D = simcase
    reffeᵤ = ReferenceFE(lagrangian, VectorValue{D,Float64}, order)
    V = TestFESpace(model, reffeᵤ, conformity=:H1,
        dirichlet_tags=["inlet", "limits", "airfoil"],
        dirichlet_masks=[(true, true), (false, true), (true, true)])
    reffeₚ = ReferenceFE(lagrangian, Float64, order)
    Q = TestFESpace(model, reffeₚ, conformity=:H1, dirichlet_tags=["outlet"])

    return V, Q
end

"""
    setup_domain_measures!(am::AirfoilModel, order::Int)

Register the volume triangulation/measure in `am.params` and return `(Ω, dΩ)`.
"""
function setup_domain_measures!(am::AirfoilModel, order::Int)
    Ω = Triangulation(am.model)
    dΩ = Measure(Ω, 2 * order)
    updatekey(am.params, :Ω, Ω)
    updatekey(am.params, :dΩ, dΩ)
    return Ω, dΩ
end

"""
    setup_shape_opt_spaces!(D, order, model, am)

Register in `am.params` the geometric quantities needed by the shape-sensitivity
kernel: interpolation space `VV0` and airfoil boundary measure/normal/tangent.
"""
function setup_shape_opt_spaces!(D::Int64, order::Int64, model, am::AirfoilModel)
    reffe = ReferenceFE(lagrangian, VectorValue{D,Float64}, order)
    VV0 = FESpace(model, reffe; conformity=:H1)

    Γ = BoundaryTriangulation(am.model; tags="airfoil")
    dΓ = Measure(Γ, 2 * order)
    nΓ = -get_normal_vector(Γ) # sign flip: Gridap normals point into the body

    updatekey(am.params, :reffe, reffe)
    updatekey(am.params, :VV0, VV0)
    updatekey(am.params, :Γ, Γ)
    updatekey(am.params, :dΓ, dΓ)
    updatekey(am.params, :nΓ, nΓ)
end

#############################################################################
# Steady solver
#############################################################################

function solve_inc_primal_steady(am::AirfoilModel, simcase::Airfoil, filename, uh00, ph00)
    @sunpack D, order, u_in, matrix_freq_update = simcase
    @unpack model = am

    V, Q = create_primal_spaces(model, simcase)

    u0 = VectorValue(u_in...)
    u_walls = VectorValue(zeros(D)...)

    U = TrialFESpace(V, [u0, u0, u_walls])
    P = TrialFESpace(Q, [0.0])

    updatekey(am.params, :U, U)
    updatekey(am.params, :P, P)

    Y = MultiFieldFESpace([V, Q])
    X = MultiFieldFESpace([U, P])

    setup_shape_opt_spaces!(D, order, model, am)
    Ω, _ = setup_domain_measures!(am, order)

    solver = LinearFESolver(LUSolver())

    uh = interpolate(u0, U)
    ph = interpolate(0.0, P)

    isnothing(uh00) || (uh.free_values .= uh00.free_values)
    isnothing(ph00) || (ph.free_values .= ph00.free_values)

    # Picard iterations: the VMS operator is re-assembled around the last iterate
    for i in 1:matrix_freq_update
        @info "Steady primal, Picard iteration $i"
        uh, ph = solve_steady_primal(uh, ph, X, Y, simcase, am.params, solver)
    end

    if !isnothing(filename)
        writevtk(Ω, filename; nsubcells=order, cellfields=["uh" => uh, "ph" => ph])
    end

    return uh, ph
end

function solve_steady_primal(uh, ph, X, Y, simcase, params, solver)
    xh = interpolate([uh, ph], X)

    updatekey(params, :uh, uh)
    res, rhs = equations_primal(simcase, params, :steady)

    op = AffineFEOperator(res, rhs, X, Y)
    Gridap.solve!(xh, solver, op)

    uh, ph = xh
    return uh, ph
end

#############################################################################
# Unsteady solver
#############################################################################

function solve_inc_primal_unsteady(am::AirfoilModel, simcase::Airfoil, filename, uh00, ph00)
    @sunpack D, order, t_endramp, t0, tF, θ, dt, u_in, time_window = simcase
    @sunpack M = simcase # save a .vtu file every M steps
    @unpack model = am

    V, Q = create_primal_spaces(model, simcase)

    u0 = ConstantInTime(VectorValue(u_in...))
    u_walls = ConstantInTime(VectorValue(zeros(D)...))
    p0 = ConstantInTime(0.0)

    U = TransientTrialFESpace(V, [u0, u0, u_walls])
    P = TransientTrialFESpace(Q, p0)

    Y = TransientMultiFieldFESpace([V, Q])
    X = TransientMultiFieldFESpace([U, P])

    setup_shape_opt_spaces!(D, order, model, am)
    Ω, _ = setup_domain_measures!(am, order)

    uh0 = interpolate(u0(0.0), U(0.0))
    ph0 = interpolate(p0(0.0), P(0.0))

    if isnothing(uh00) && t_endramp == t0
        # initialize from the steady solution
        uh00, ph00 = solve_inc_primal_steady(am, simcase, nothing, uh00, uh00)
    end

    uh0.free_values .= uh00.free_values
    ph0.free_values .= ph00.free_values

    xh0 = interpolate([uh0, ph0], X(0.0))

    updatekey(am.params, :uh, uh0)
    updatekey(am.params, :ph, ph0)

    m, res, rhs = equations_primal(simcase, am.params, :unsteady)
    op = TransientLinearFEOperator((res, m), rhs, X, Y)

    ode_solver = ThetaMethod(LUSolver(), dt, θ)
    sol = Gridap.solve(ode_solver, op, t0, tF, xh0)

    UH = [copy(uh0.free_values)]
    PH = [copy(ph0.free_values)]

    res_path = "Results_unsteady_primal"
    mkpath(res_path)

    createpvd(filename) do pvd
        pvd[0] = createvtk(Ω, joinpath(res_path, "$(filename)_0.vtu");
            nsubcells=order, cellfields=["uh" => uh0, "ph" => ph0])
        for (idx, (t, xhtn)) in enumerate(sol)
            uh, ph = xhtn
            push!(UH, copy(uh.free_values))
            push!(PH, copy(ph.free_values))
            println("Primal solved at time step $t")

            # keep the linearization field of the VMS operator at the current state
            copyto!(am.params[:uh].free_values, uh.free_values)

            if mod(idx, M) == 0
                pvd[t] = createvtk(Ω, joinpath(res_path, "$(filename)_$t.vtu");
                    nsubcells=order, cellfields=["uh" => uh, "ph" => ph])
            end
        end
    end

    updatekey(am.params, :UH, UH)
    updatekey(am.params, :PH, PH)

    jldsave("UnsteadyPrimalFields.jld2"; am)

    avg_UH, avg_PH = time_average_fields(UH, PH, time_window, dt, t0)
    uh0.free_values .= avg_UH
    ph0.free_values .= avg_PH

    return uh0, ph0
end

"""
    time_average_fields(UH, PH, time_window, dt, t0)

Average the snapshot vectors `UH`, `PH` over the steps whose time lies in `time_window`.
"""
function time_average_fields(UH, PH, time_window, dt::Float64, t0::Float64)
    t_lo, t_hi = time_window
    times = t0 .+ dt .* (0:length(UH)-1)
    indices = findall(t -> t_lo ≤ t ≤ t_hi, times)

    return Statistics.mean(UH[indices]), Statistics.mean(PH[indices])
end
