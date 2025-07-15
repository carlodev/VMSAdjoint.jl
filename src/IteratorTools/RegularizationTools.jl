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



function generate_regularized_model(adesign::AirfoilDesign, i::Int64, ss::Float64, meshinfo, physicalp, folder::String; initial_R=0.0, max_tries=50)
    i_try = 0
    model = nothing
    R = initial_R
    flag = true

    function regf(x0, y0)
        y1 = y0
        R > 0.0 && println("Denoise Radius $R")
        R > 0 && (y1, _ = denoise(y0; factor=R))
        return y1
    end

    reg = Regularization(active=true, iter_reg=1, fun=regf)

    while flag && i_try < max_tries
        adesign_tmp = adesign
        if ss> 0.0 
            adesign_tmp = perturb_DesignParameter(adesign, i, ss)
        end

        adesign_r = regularize_airfoil(adesign_tmp, 1, reg)
        modelname = create_msh(meshinfo, adesign_r, physicalp, folder; iter=i)
        
        try
            model = GmshDiscreteModel(modelname)
        catch
            i_try += 1
            R += 0.01
            println("Mesh gen $(i_try)")
        else
            flag = false
        end
    end

    return model
end