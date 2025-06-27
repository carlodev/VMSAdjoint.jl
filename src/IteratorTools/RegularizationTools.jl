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