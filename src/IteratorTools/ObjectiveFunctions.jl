#############################################################################
# Aerodynamic coefficients, objective function and its derivative
#############################################################################

"""
    compute_airfoil_forces(uh, ph, nΓ, dΓ, ν)

Drag and lift over the airfoil boundary, accounting for pressure and viscous
stresses. `nΓ` must point outward with respect to the body.
"""
function compute_airfoil_forces(uh::SingleFieldFEFunction, ph::SingleFieldFEFunction, nΓ::OperationCellField, dΓ::GenericMeasure, ν::Float64)
    IForce = ∫(-ph ⋅ nΓ + ν * (∇(uh) + transpose(∇(uh))) ⋅ nΓ)dΓ
    D, L = sum(IForce)
    return D, L
end

"""
    compute_airfoil_coefficients(uh, ph, nΓ, dΓ, physicalp)

Normalize the airfoil forces by the dynamic-pressure reference, obtaining `(CD, CL)`.
"""
function compute_airfoil_coefficients(uh::SingleFieldFEFunction, ph::SingleFieldFEFunction, nΓ::OperationCellField, dΓ::GenericMeasure, physicalp::PhysicalParameters)
    @unpack c, u_in_mag, ν = physicalp

    q = 0.5 * c * u_in_mag^2 # dynamic pressure reference C∞ (ρ=1); must match compute_gradient

    D, L = compute_airfoil_forces(uh, ph, nΓ, dΓ, ν)
    CD = D / q
    CL = L / q

    @info "CL = $CL ; CD = $CD"

    return CD, CL
end

"""
    obj_fun(am::AirfoilModel, vbcase::Airfoil, uh, ph, thick_penalty, fun::Function)

Value of the objective function `fun([CD, CL])`, augmented with the
minimum-thickness penalty. Returns `(fitness + penalty, [CD, CL])`.
"""
function obj_fun(am::AirfoilModel, vbcase::Airfoil, uh, ph, thick_penalty::ThickPenalty, fun::Function)
    @sunpack order = vbcase
    @unpack model = am
    @unpack thickness_penalty, valid = thick_penalty

    Γ = BoundaryTriangulation(model; tags="airfoil")
    dΓ = Measure(Γ, 2 * order)
    nΓ = -1 .* get_normal_vector(Γ)

    physicalp = vbcase.simulationp.physicalp
    CD, CL = compute_airfoil_coefficients(uh, ph, nΓ, dΓ, physicalp)

    thick_pen = valid ? thickness_penalty(interpolate_points_x(am)) : 0.0
    fitnessval = fun([CD, CL])

    @info "fitnessval = $fitnessval ; thickness_penalty = $thick_pen"

    return fitnessval + thick_pen, [CD, CL]
end

"""
    dJobj_fun(fun::Function, CDCL::Vector{Float64})

Gradient of the objective function with respect to `[CD, CL]` via automatic
differentiation; its negative is the adjoint boundary condition on the airfoil.
"""
function dJobj_fun(fun::Function, CDCL::Vector{Float64})
    @assert length(CDCL) == 2
    ForwardDiff.gradient(fun, CDCL)
end

#############################################################################
# Thickness sampling
#############################################################################

"Filter one airfoil surface to `x ∈ [xx0, xx1]`, asserting sortedness."
function _filter_surface(x::Vector{Float64}, y::Vector{Float64}, xx0::Float64, xx1::Float64, side::String)
    valid_idx = findall(xi -> xx0 <= xi <= xx1, x)
    x_f, y_f = x[valid_idx], y[valid_idx]
    @assert issorted(x_f) "$side surface x-coordinates not sorted"
    return x_f, y_f
end

"""
    interpolate_points_x(am::AirfoilModel)

Sample both airfoil surfaces on a common x-grid (excluding small leading/trailing
edge cutoffs), for the thickness-penalty evaluation. Returns `(xx, yu, yl)`.
"""
function interpolate_points_x(am::AirfoilModel)
    leading_edge_cutoff = 0.01
    trailing_edge_cutoff = 0.01

    # overall valid x-range of the airfoil
    xx0 = maximum([minimum(am.ap.xu); minimum(am.ap.xl)]) + leading_edge_cutoff
    xx1 = minimum([maximum(am.ap.xu); maximum(am.ap.xl)]) - trailing_edge_cutoff
    @assert xx1 > xx0 "Airfoil x-coordinates don't span a valid range"

    xu_f, yu_f = _filter_surface(am.ap.xu, am.ap.yu, xx0, xx1, "Upper")
    xl_f, yl_f = _filter_surface(am.ap.xl, am.ap.yl, xx0, xx1, "Lower")

    # common evaluation points
    xx01 = maximum([minimum(xu_f); minimum(xl_f)])
    xx11 = minimum([maximum(xu_f); maximum(xl_f)])
    xx = collect(LinRange(xx01, xx11, 201))

    yu = linear_interpolation(xu_f, yu_f).(xx)
    yl = linear_interpolation(xl_f, yl_f).(xx)

    return (xx, yu, yl)
end
