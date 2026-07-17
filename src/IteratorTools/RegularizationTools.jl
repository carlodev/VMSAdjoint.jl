function regularize_airfoil(ad::AirfoilDesign, iter::Int64, regularization::Regularization)
    @unpack active, iter_reg, fun = regularization

    ad_reg = active && mod(iter, iter_reg) == 0 ? regularize_airfoil(ad, fun) : ad

    return ad_reg
end

function regularize_airfoil(ad::AirfoilDesign, fun::Function)
    println("--- Regularization Airfoil ---")
    yu_reg = fun(ad.ap.xu, ad.ap.yu) #regularize top points airfoil
    yl_reg = fun(ad.ap.xl, ad.ap.yl) #regularize bottom points airfoil
    ap_new = AirfoilPoints(ad.ap.xu, ad.ap.xl, yu_reg, yl_reg)
    regularize_airfoil(ad, ap_new)
end

function regularize_airfoil(ad::RBFDesign, ap_new::AirfoilPoints)
    rbfd_new = RBFDesign(ad.rbfg, ap_new)
    return rbfd_new
end


function regularize_airfoil(ad::AirfoilCSTDesign, ap_new::AirfoilPoints)
    @error "Regularization when using CST not advised !"
end



"""
    generate_model(adesign, i, ss, meshinfo, physicalp, folder; max_cells=1_500_000)

Build a `GmshDiscreteModel` from `adesign`, perturbing design parameter `i` by the
signed shift `ss` (`ss = 0` -> baseline geometry). No surface regularization is
applied here: design regularization happens only when the user explicitly requests
it through `AdjSolver(regularization = ...)`, which is consumed by
`regularize_airfoil` before this function is called.

`max_cells` guards against runaway mesh density (seen with mismatched Gmsh
versions, where background-field sizing behaves differently): the run aborts
immediately instead of hanging in the linear solver.
"""
function generate_model(adesign::AirfoilDesign, i::Int64, ss::Float64, meshinfo, physicalp, folder::String; max_cells::Int=1_500_000)
    adesign_tmp = ss != 0.0 ? perturb_DesignParameter(adesign, i, ss) : adesign

    modelname = create_msh(meshinfo, adesign_tmp, physicalp, folder; iter=i)

    model = GmshDiscreteModel(modelname)

    ncells = num_cells(model)
    ncells <= max_cells || error(
        "Generated mesh has $ncells cells (> max_cells = $max_cells). " *
        "This usually means the Gmsh version differs from the one the sizing was tuned on " *
        "(check `gmsh.option.getString(\"General.Version\")`). " *
        "Raise `max_cells` only if the density is intentional.")

    return model
end