function perturb_DesignParameter(cstd::AirfoilCSTDesign, i::Int64, ss::Float64 )
    Ndes = length(cstd.cstg.cstw)
    @assert i <= Ndes
    N = length(cstd.ap.xu) + length(cstd.ap.xl) - 1
    cstg_new = perturb_CSTweights(cstd.cstg, i, ss)
    cstd_new = AirfoilCSTDesign(cstg_new,N) 
    return cstd_new
end

function perturb_CSTweights(cstg::CSTGeometry, i::Int, ss::Real)
    cstw_new = perturb_CSTweights(cstg.cstw, i, ss) 
    cstg_new = CSTGeometry(cstw_new, cstg.dz,cstg.N1,cstg.N2)
    return cstg_new
end

function perturb_CSTweights(cstw::CSTweights, i::Int, ss::Real)
    nu = length(cstw.wu)

    new_wu = copy(cstw.wu)
    new_wl = copy(cstw.wl)

    if i <= nu
        new_wu[i] += ss
    else
        new_wl[i - nu] += ss
    end

    return CSTweights(new_wu, new_wl)
end


function AirfoilTools.AirfoilCST.CSTweights(w::Vector{Float64})
    Ndes = length(w)
    @assert iseven(Ndes) "The number of design parameters must be even"
    Nhalf = Int(Ndes ÷ 2)
    wu = w[1:Nhalf]
    wl = w[Nhalf+1:end]
    return CSTweights(wu,wl)
end

function AirfoilCSTDesign(cstd::AirfoilCSTDesign,w::Vector{Float64})

    N = length(cstd.ap.xu) + length(cstd.ap.xl) - 1

    cstw_new = CSTweights(w)
    cstg_new = CSTGeometry(cstw_new, cstd.cstg.dz,cstd.cstg.N1,cstd.cstg.N2)
    cstd_new = AirfoilCSTDesign(cstg_new,N) 
    
    return cstd_new
end
