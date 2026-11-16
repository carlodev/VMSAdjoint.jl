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
# Dimension-generic boundary-condition helpers
#############################################################################

"Dirichlet mask constraining every velocity component."
full_mask(D::Int) = ntuple(_ -> true, D)

"Dirichlet mask constraining only the wall-normal (y) component: free-slip lateral boundaries."
wall_normal_mask(D::Int) = ntuple(i -> i == 2, D)

"""
    boundary_vector(vals, D)

`VectorValue{D}` from the (possibly shorter) vector `vals`, zero-padding the
remaining components. Used to build D-dimensional Dirichlet values from the
2-component adjoint boundary condition `-dJ/d[CD,CL]` (the spanwise component
is zero in 3D).
"""
boundary_vector(vals::AbstractVector{<:Real}, D::Int) =
    VectorValue(ntuple(i -> i <= length(vals) ? Float64(vals[i]) : 0.0, D))

"""
    set_free_values!(fh, fh0, name)

Copy the free values of `fh0` into `fh` when `fh0` is provided and the DOF
counts match; warn and keep `fh` unchanged on a mismatch (e.g. warm-starting
across a remeshing that changed the number of DOFs). Returns whether the copy
happened.
"""
function set_free_values!(fh, fh0, name::String)
    isnothing(fh0) && return false
    n, n0 = length(fh.free_values), length(fh0.free_values)
    if n == n0
        fh.free_values .= fh0.free_values
        return true
    end
    @warn "Skipping initial values for $name: DOF count changed ($n0 -> $n)"
    return false
end

#############################################################################
# Entry points
#############################################################################

"""
    solve_inc_primal(am, simcase, [filename,] timed::Symbol; uh0=nothing, ph0=nothing, linear_backend=:cpu, picard_iters=nothing)

Solve the primal incompressible flow problem, `timed ∈ (:steady, :unsteady)`.
When `filename` is omitted it defaults to `"inc-primal-<timed>"`.
`linear_backend ∈ (:cpu, :gpu)` selects the linear solver (see [`linear_solver`](@ref)).
Returns `(uh, ph)`; for `:unsteady` these are the time-averaged fields over
`time_window`, and the full history is stored in `am.params[:UH]/[:PH]`.

`picard_iters` overrides how many Picard sweeps the STEADY solver takes. It
exists because the two solvers want opposite things from
`SolverParameters.matrix_freq_update`, which used to drive both: the transient
solver wants it SMALL (a stale linearization is refreshed often), the steady
solver wants it LARGE (that is its iteration count). Any run that needs a
converged steady state *and* a cheap transient — E4's steady/steady leg does —
must set the two independently. `nothing` keeps the historical behaviour of
taking the count from `matrix_freq_update`; ignored when `timed = :unsteady`.
"""
function solve_inc_primal(am::AirfoilModel, simcase::Airfoil, timed::Symbol; uh0=nothing, ph0=nothing, linear_backend::Symbol=:cpu, picard_iters::Union{Nothing,Int}=nothing)
    solve_inc_primal(am, simcase, "inc-primal-$timed", timed; uh0=uh0, ph0=ph0, linear_backend=linear_backend, picard_iters=picard_iters)
end

function solve_inc_primal(am::AirfoilModel, simcase::Airfoil, filename::String, timed::Symbol; uh0=nothing, ph0=nothing, linear_backend::Symbol=:cpu, picard_iters::Union{Nothing,Int}=nothing)
    if timed === :steady
        solve_inc_primal_steady(am, simcase, filename, uh0, ph0, linear_backend, picard_iters)
    elseif timed === :unsteady
        solve_inc_primal_unsteady(am, simcase, filename, uh0, ph0, linear_backend)
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
        dirichlet_masks=[full_mask(D), wall_normal_mask(D), full_mask(D)])
    reffeₚ = ReferenceFE(lagrangian, Float64, order)
    Q = TestFESpace(model, reffeₚ, conformity=:H1, dirichlet_tags=["outlet"])

    return V, Q
end

"""
    velocity_pressure_space(Vu, Qp, backend)

Product velocity/pressure space, block-assembled for the `:block` backend and
interleaved (Gridap's default) otherwise — see [`multifield_style`](@ref).
Works for both plain and transient component spaces, since
`TransientMultiFieldFESpace === MultiFieldFESpace`.
"""
velocity_pressure_space(Vu, Qp, backend::Symbol) =
    MultiFieldFESpace([Vu, Qp]; style=multifield_style(backend))

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
kernel: interpolation space `VV0` and airfoil boundary measure/normal (via
[`Interfaces.airfoil_boundary!`](@ref), shared with the force computation).
"""
function setup_shape_opt_spaces!(D::Int64, order::Int64, model, am::AirfoilModel)
    reffe = ReferenceFE(lagrangian, VectorValue{D,Float64}, order)
    VV0 = FESpace(model, reffe; conformity=:H1)

    updatekey(am.params, :reffe, reffe)
    updatekey(am.params, :VV0, VV0)
    airfoil_boundary!(am, order)
end

#############################################################################
# Steady solver
#############################################################################

function solve_inc_primal_steady(am::AirfoilModel, simcase::Airfoil, filename, uh00, ph00, linear_backend::Symbol=:cpu, picard_iters::Union{Nothing,Int}=nothing)
    @sunpack D, order, u_in, matrix_freq_update = simcase
    npicard = isnothing(picard_iters) ? matrix_freq_update : picard_iters
    npicard >= 1 || throw(ArgumentError("picard_iters must be >= 1, got $npicard"))
    @unpack model = am

    V, Q = create_primal_spaces(model, simcase)

    u0 = VectorValue(u_in...)
    u_walls = VectorValue(zeros(D)...)

    U = TrialFESpace(V, [u0, u0, u_walls])
    P = TrialFESpace(Q, [0.0])

    updatekey(am.params, :U, U)
    updatekey(am.params, :P, P)

    Y = velocity_pressure_space(V, Q, linear_backend)
    X = velocity_pressure_space(U, P, linear_backend)

    setup_shape_opt_spaces!(D, order, model, am)
    Ω, _ = setup_domain_measures!(am, order)

    # built after the measures: the :block backend assembles a pressure Schur
    # approximation over dΩ
    ctx = solver_context(simcase, am.params, P, Q; transient=false)
    solver = LinearFESolver(linear_solver(linear_backend, ctx))

    uh = interpolate(u0, U)
    ph = interpolate(0.0, P)

    set_free_values!(uh, uh00, "uh (steady init)")
    set_free_values!(ph, ph00, "ph (steady init)")

    # Picard iterations: the VMS operator is re-assembled around the last iterate.
    # `npicard` is a fixed iteration count, not a convergence tolerance, so the
    # relative increment is logged every sweep to make silent under-convergence
    # visible instead of assuming N iterations was enough. The loop stops early
    # once the tolerance is met — with `picard_iters` set generously that is the
    # normal exit, and the sweeps after convergence would be pure waste.
    picard_tol = 1e-6
    rel_du, rel_dp = Inf, Inf
    for i in 1:npicard
        uh_prev, ph_prev = uh, ph
        @time uh, ph = solve_steady_primal(uh, ph, X, Y, simcase, am.params, solver)

        rel_du = norm(get_free_dof_values(uh) .- get_free_dof_values(uh_prev)) / norm(get_free_dof_values(uh))
        rel_dp = norm(get_free_dof_values(ph) .- get_free_dof_values(ph_prev)) / norm(get_free_dof_values(ph))
        @info "Steady primal, Picard iteration $i/$npicard: rel Δu = $rel_du, rel Δp = $rel_dp"

        max(rel_du, rel_dp) <= picard_tol && break
    end

    if max(rel_du, rel_dp) > picard_tol
        @warn "Steady primal did not reach the Picard convergence tolerance ($picard_tol) " *
              "after $npicard iterations: rel Δu = $rel_du, rel Δp = $rel_dp. " *
              "Raise `picard_iters` if this run feeds a gradient comparison; if the residual " *
              "has stalled rather than decayed, the steady state is likely unstable at this " *
              "operating point and Picard will not reach it however many sweeps it is given."
    end
    updatekey(am.params, :picard_residual, (rel_du=rel_du, rel_dp=rel_dp, tol=picard_tol))

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

function solve_inc_primal_unsteady(am::AirfoilModel, simcase::Airfoil, filename, uh00, ph00, linear_backend::Symbol=:cpu)
    @sunpack D, order, t0, tF, θ, dt, u_in, time_window = simcase
    @sunpack M = simcase # save a .vtu file every M steps
    @sunpack matrix_freq_update = simcase # rebuild+refactorize the operator every N steps
    @unpack model = am

    V, Q = create_primal_spaces(model, simcase)

    u0 = ConstantInTime(VectorValue(u_in...))
    u_walls = ConstantInTime(VectorValue(zeros(D)...))
    p0 = ConstantInTime(0.0)

    U = TransientTrialFESpace(V, [u0, u0, u_walls])
    P = TransientTrialFESpace(Q, p0)

    Y = velocity_pressure_space(V, Q, linear_backend)
    X = velocity_pressure_space(U, P, linear_backend)

    setup_shape_opt_spaces!(D, order, model, am)
    Ω, _ = setup_domain_measures!(am, order)

    uh0 = interpolate(u0(0.0), U(0.0))
    ph0 = interpolate(p0(0.0), P(0.0))

    # Impulsive start from the freestream interpolant: the start-up transient is
    # what seeds the wake instability. Initializing from the converged steady
    # solution (a Picard solve) instead would place the flow exactly on the
    # unstable equilibrium with no perturbation to grow, pinning it at steady
    # state and suppressing vortex shedding. A caller may still warm-start from a
    # supplied field by passing uh00/ph00.
    set_free_values!(uh0, uh00, "uh (unsteady init)")
    set_free_values!(ph0, ph00, "ph (unsteady init)")

    xh0 = interpolate([uh0, ph0], X(t0))

    updatekey(am.params, :uh, uh0)
    updatekey(am.params, :ph, ph0)

    UH = [copy(uh0.free_values)]
    PH = [copy(ph0.free_values)]

    res_path = "Results_unsteady_primal"
    mkpath(res_path)

    ctx = solver_context(simcase, am.params, P, Q; transient=true)
    lsolver = linear_solver(linear_backend, ctx)

    # rebuilt once per chunk, around the linearization field currently in
    # params[:uh]; `constant_forms` is what lets Gridap reuse the factorization
    # for every step inside the chunk
    make_operator = function ()
        m, res, rhs = equations_primal(simcase, am.params, :unsteady)
        return TransientLinearFEOperator((res, m), rhs, X, Y; constant_forms=(true, true))
    end

    # advance the frozen linearization to the latest velocity state
    refresh! = function (_step, state)
        xh_state = FEFunction(X(t0), state)
        copyto!(am.params[:uh].free_values, get_free_dof_values(first(xh_state)))
    end

    createpvd(filename) do pvd
        pvd[t0] = createvtk(Ω, joinpath(res_path, "$(filename)_0.vtu");
            nsubcells=order, cellfields=["uh" => uh0, "ph" => ph0])

        on_step! = function (idx, t, xhtn)
            uh, ph = xhtn
            push!(UH, copy(uh.free_values))
            push!(PH, copy(ph.free_values))

            if mod(idx, M) == 0
                pvd[t] = createvtk(Ω, joinpath(res_path, "$(filename)_$t.vtu");
                    nsubcells=order, cellfields=["uh" => uh, "ph" => ph])
            end
        end

        march_frozen_theta(make_operator, lsolver, X, xh0, t0, tF, dt, θ, matrix_freq_update;
            on_step! = on_step!, refresh! = refresh!, label = "Primal")
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
