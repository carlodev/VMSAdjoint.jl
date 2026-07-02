
"""
    compute_sensitivity(am0, am1, adesign, IIDX, δ, simcase, thick_penalty, uh, uhadj)

From the primal solution `uh` and the adjoint solution `uhadj`, compute the shape
sensitivity of the (total-force) objective with respect to design parameter `IIDX`.

am0 :: AirfoilModel of the un-perturbed geometry
am1 :: AirfoilModel of the perturbed geometry
δ   :: signed perturbation of the design parameter

The two supported parametrizations (RBF and CST) differ *only* in how the normal
mesh-deformation field is built; that part is dispatched through
`deformation_normal_field`, while the actual gradient kernel lives in
`compute_gradient`.
"""
function compute_sensitivity(am0::AirfoilModel, am1::AirfoilModel, adesign::AirfoilDesign, IIDX::Int64, δ::Float64, simcase::Airfoil, thick_penalty::ThickPenalty, uh, uhadj)
    @sunpack ν = simcase
    @unpack tΓ, nΓ, dΓ = am0.params

    physicalp = simcase.simulationp.physicalp
    C∞ = 0.5 * physicalp.c * physicalp.u_in_mag^2  # force-coefficient reference; must match ObjectiveFunctions.q

    v_field_corr = deformation_normal_field(adesign, IIDX, am0, am1, δ, simcase)

    return compute_gradient(uh, uhadj, ν, tΓ, nΓ, dΓ, am0, am1, δ, thick_penalty, v_field_corr, C∞)
end


"""
    deformation_normal_field(adesign, IIDX, am0, am1, δ, simcase)

Normal component (δβ) on the airfoil boundary of the mesh-deformation velocity
induced by perturbing design parameter `IIDX`. One method per parametrization.
"""
function deformation_normal_field(::AirfoilCSTDesign, IIDX::Int64, am0::AirfoilModel, am1::AirfoilModel, δ::Float64, simcase::Airfoil)
    @unpack VV0, nΓ, reffe = am0.params

    fx(x) = x
    m0 = get_free_dof_values(interpolate_everywhere(fx, VV0))          # node positions, un-perturbed mesh
    VV1 = FESpace(am1.model, reffe; conformity=:H1)
    m1 = get_free_dof_values(interpolate_everywhere(fx, VV1))          # node positions, perturbed mesh
    vi = (m1 - m0) ./ δ                                                # mesh velocity ≈ ∂x/∂β

    v_field = FEFunction(VV0, vi)
    return v_field ⋅ nΓ
end


function deformation_normal_field(rbfd::RBFDesign, IIDX::Int64, am0::AirfoilModel, am1::AirfoilModel, δ::Float64, simcase::Airfoil)
    @sunpack AoA = simcase
    @unpack VV0, nΓ = am0.params

    function rotation(a::Vector{Float64}, AoA::Float64)
        x, y = a
        xrt = x * cosd(AoA) + y * sind(AoA)
        yrt = -1 * x * sind(AoA) + y * cosd(AoA)
        return [xrt, yrt]
    end

    top_cp(n::VectorValue{2,Float64})    = n[2] > 0 ? VectorValue(0.0,  1.0) : VectorValue(0.0, 0.0)
    bottom_cp(n::VectorValue{2,Float64}) = n[2] < 0 ? VectorValue(0.0, -1.0) : VectorValue(0.0, 0.0)

    Ndes = length(get_DesignParameters(rbfd))
    Nhalf = Int(Ndes ÷ 2)

    xr, yr = rotation([rbfd.rbfg.control_points[IIDX], rbfd.cy[IIDX]], AoA)
    RB = rbfd.rbfg.RBFfun

    function fx_rbf(x)
        x1, y1 = rotation([x...], -AoA)
        dist = norm([x1, y1] - [xr, yr])
        rb = AirfoilTools.AirfoilRBF.fRBF(dist, RB.support_radius, RB.fun)
        return VectorValue(0.0, float(rb))
    end

    vv_corr = IIDX <= Nhalf ? top_cp ∘ nΓ : bottom_cp ∘ nΓ

    v_field = interpolate_everywhere(fx_rbf, VV0)
    return v_field ⋅ vv_corr
end


"""
    compute_gradient(uh, uhadj, ν, tΓ, nΓ, dΓ, am0, am1, δ, thick_penalty, v_field_corr, C∞)

Higher-derivative-free surface sensitivity (Sorgiovanni 2016, Eq. 3.78 / Castro et al.,
Total-Force LFA):

    dJ/dβ = (1/C∞) ∫_Γ  ν S δβ  dΓ ,   S = ⟨(∂v_t/∂n)(∂ψ_t/∂n)⟩

with v = primal velocity (uh), ψ = adjoint velocity (uhadj), t/n the wall tangent/normal
and δβ = `v_field_corr` the normal mesh displacement. In Gridap's convention
(∇u)[i,j] = ∂u_j/∂x_i, the contraction `(∇(u)⋅tΓ)⋅nΓ` evaluates to ∂u_t/∂n.

- Steady / averaged flow: S is the instantaneous product (∂v_t/∂n)(∂ψ_t/∂n).
- Unsteady flow: the unsteady adjoint solver builds the *time-averaged product*
  ⟨(∂v_t/∂n)(∂ψ_t/∂n)⟩ (which keeps the primal–adjoint covariance term that a
  product of time-averages would drop — see Srinath & Mittal, JCP 2010, Eq. 12) and
  stores it in `am0.params[:wall_shear_corr]`. When present, it is used directly.

`C∞` is passed explicitly (previously hard-coded as the factor 2, valid only for
c = u∞ = 1). The kernel returns the same value as before when C∞ = 0.5.
"""
function compute_gradient(uh, uhadj, ν, tΓ, nΓ, dΓ, am0, am1, δ, thick_penalty, v_field_corr, C∞)
    ∂ₜ∂n(u) = (∇(u) ⋅ tΓ) ⋅ nΓ  # ∂u_t/∂n

    # unsteady: time-averaged wall-shear correlation built by the adjoint solver;
    # steady: instantaneous product of the two supplied fields.
    S = get(am0.params, :wall_shear_corr, nothing)
    wall_shear = isnothing(S) ? ∂ₜ∂n(uh) * ∂ₜ∂n(uhadj) : S

    J_sens = (1 / C∞) * sum(∫(ν * wall_shear * v_field_corr)dΓ)

    # thickness-constraint penalty gradient (finite difference between the two designs)
    Jp = compute_∇penalty(am0, am1, δ, thick_penalty)

    return J_sens, Jp
end


function compute_∇penalty(am0::AirfoilModel, am1::AirfoilModel, δ::Float64, thick_penalty::ThickPenalty)
    @unpack thickness_penalty, valid = thick_penalty

    interpolated_points0 = interpolate_points_x(am0)
    interpolated_points1 = interpolate_points_x(am1)

    penalty0 = valid ? thickness_penalty(interpolated_points0) : 0.0
    penalty1 = valid ? thickness_penalty(interpolated_points1) : 0.0

    ∇penalty = (penalty1 - penalty0) / δ

    return ∇penalty
end
