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



function generate_regularized_model(adesign::AirfoilDesign, i::Int64, ss::Float64, meshinfo, physicalp, folder::String; initial_L=0.0, max_tries=10)
    i_try = 0
    model = nothing
    L = initial_L
    flag = true

    function regf(x0, y0)
        y1 = deepcopy(y0)
        L > 0.0 && println("Denoise L $L")
        # R > 0 && (y1, _ = denoise(y0; factor=R))
        fn = fit_sine_series(x0,y0, 25,lambda = L) 
        L > 0 && (y1 = fn.(x0) ) 
        y1[1:2] = y0[1:2]
        y1[end-1:end] = y0[end-1:end]

        return y1
    end

    reg = Regularization(active=true, iter_reg=1, fun=regf)

    while flag && i_try <= max_tries
        adesign_tmp = adesign
        if ss> 0.0 
            adesign_tmp = perturb_DesignParameter(adesign, i, ss)
        end

        adesign_r = regularize_airfoil(adesign_tmp, 1, reg)
        modelname = create_msh(meshinfo, adesign_r, physicalp, folder; iter=i)
        
        try
            model = GmshDiscreteModel(modelname)
            flag = false
        catch
            i_try += 1
            L += 0.00001
            println("Mesh gen $(i_try)")
            if i_try==max_tries
                @error "Impossible to generate a regularized model"
            end

        end
    end

    return model
end