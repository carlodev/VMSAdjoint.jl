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



# Airfoil-surface denoising closure used to smooth a (possibly perturbed) design
# before meshing. `L` is the sine-series regularization weight; TE/LE points are clamped.
_denoise_fun(L::Float64) = function (x0, y0)
    L <= 0.0 && return deepcopy(y0)
    println("Denoise L $L")
    y1 = fit_sine_series(x0, y0, 25; lambda=L).(x0)
    y1[1:2] .= y0[1:2]
    y1[end-1:end] .= y0[end-1:end]
    return y1
end

"""
    generate_regularized_model(adesign, i, ss, meshinfo, physicalp, folder; initial_L=0.0, max_tries=10, L_step=1e-5)

Build a `GmshDiscreteModel` from `adesign`. If Gmsh produces a mesh that GridapGmsh
cannot read, progressively increase the surface-smoothing weight `L` and retry, up to
`max_tries` times. Gmsh is finalized before each attempt so a failed run never leaks
state into the next one. Only mesh-generation failures are retried; interrupts and
out-of-memory errors are rethrown.
"""
function generate_regularized_model(adesign::AirfoilDesign, i::Int64, ss::Float64, meshinfo, physicalp, folder::String; initial_L=0.0, max_tries=10, L_step=1e-5)
    L = initial_L

    for i_try in 0:max_tries
        # tear down any Gmsh session left open by a previous failed attempt
        gmsh.isInitialized() == 1 && gmsh.finalize()

        adesign_tmp = ss > 0.0 ? perturb_DesignParameter(adesign, i, ss) : adesign
        reg = Regularization(active=true, iter_reg=1, fun=_denoise_fun(L))
        adesign_r = regularize_airfoil(adesign_tmp, 1, reg)
        modelname = create_msh(meshinfo, adesign_r, physicalp, folder; iter=i)

        try
            return GmshDiscreteModel(modelname)
        catch err
            err isa InterruptException && rethrow()
            @warn "Mesh generation failed (try $i_try), increasing smoothing" L exception = err
            L += L_step
        end
    end

    error("Impossible to generate a regularized model after $max_tries tries")
end