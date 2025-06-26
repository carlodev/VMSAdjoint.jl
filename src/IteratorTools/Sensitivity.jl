
"""
    compute_sensitivity(am0,am1,  δ::Float64, simcase::Airfoil,  uh,ph,uhadj,phadj)

From the solution of the primal flow uh,ph and the adjoint solution uhadj and phadj it computes the senstivities.
am0:: Airfoil model of un-pertubed geometry
am0:: Airfoil model of pertubed geometry
δ:: is the signed perturbation
"""

function compute_sensitivity(am0::AirfoilModel,am1::AirfoilModel, adesign::AirfoilCSTDesign, IIDX::Int64, δ::Float64, simcase::Airfoil,thick_penalty::ThickPenalty,  uh,uhadj)

    @sunpack  ν = simcase
    @unpack VV0, dΩ, reffe = am0.params
    @unpack tΓ, nΓ, dΓ = am0.params

    function fx(x)
        return x
    end
    
    m0 = get_free_dof_values(interpolate_everywhere(fx,VV0)) #it gives the nodes position in the un-perturbed mesh
    VV1 = FESpace(am1.model, reffe; conformity=:H1)
    m1 = get_free_dof_values(interpolate_everywhere(fx,VV1)) #it gives the nodes position in the perturbed mesh
    vi = (m1-m0)./ δ

    v_field = FEFunction(VV0,vi)
   

    v_field_corr = v_field ⋅nΓ

    compute_gradient(uh, uhadj, ν, tΓ, nΓ,dΓ, am0,am1,  δ, thick_penalty, v_field_corr )
    
end






function compute_sensitivity(am0::AirfoilModel,am1::AirfoilModel, rbfd::RBFDesign, IIDX::Int64, δ::Float64, simcase::Airfoil,thick_penalty::ThickPenalty,  uh,uhadj)

    @sunpack ν, AoA = simcase
    @unpack VV0, dΩ, reffe = am0.params
    @unpack tΓ, nΓ, dΓ = am0.params

    function rotation(a::Vector{Float64},AoA::Float64)
        x,y = a
        xrt= x * cosd(AoA) + y*sind(AoA)
        yrt = -1*x * sind(AoA) + y *cosd(AoA)
        return [xrt,yrt]
    end

    function top_cp(n::VectorValue{2,Float64})
        n1, n2 = n
        n2 > 0 ? VectorValue(0.0, 1.0) : VectorValue(0.0, 0.0) 
    end
    
    function bottom_cp(n::VectorValue{2,Float64})
        n1, n2 = n
        n2 < 0 ? VectorValue(0.0, -1.0) : VectorValue(0.0, 0.0) 
    end
    
    

    Ndes = length(get_DesignParameters(rbfd))
    Nhalf = Int(Ndes ÷ 2)
    
    xr,yr=rotation([rbfd.rbfg.control_points[IIDX],rbfd.cy[IIDX]],AoA)
    
    RB = rbfd.rbfg.RBFfun

    function fx_rbf(x)
            x1,y1 = rotation([x...],-AoA)
            if 0.0<=x1<0.98
                dist = norm([x1,y1] - [xr,yr] )
                rb = AirfoilTools.AirfoilRBF.fRBF(dist, RB.support_radius, RB.fun)
            else
                rb = 0.0
            end

            return  VectorValue(0.0,float(rb))
    end

        if IIDX <= Nhalf
            vv_corr = top_cp ∘ nΓ
        else
            vv_corr = bottom_cp  ∘ nΓ
        end



    v_field = interpolate_everywhere(fx_rbf,VV0)
    v_field_corr = v_field⋅vv_corr
    compute_gradient(uh, uhadj, ν, tΓ, nΓ,dΓ, am0,am1,  δ, thick_penalty, v_field_corr )
    
end


function compute_gradient(uh, uhadj, ν, tΓ, nΓ,dΓ, am0,am1,  δ, thick_penalty, v_field_corr )
    J_sens = -2 .*    sum(-ν* ∫( ((∇(uh) ⋅ tΓ) ⋅ nΓ) * ((∇(uhadj) ⋅ tΓ) ⋅ nΓ) ⋅ v_field_corr)dΓ)
    #ThickPenalty Gradient
    Jp = compute_∇thickness_penalty(am0,am1,  δ, thick_penalty)
return  J_sens,Jp
end






function compute_∇thickness_penalty(am0::AirfoilModel,am1::AirfoilModel,  δ::Float64, thick_penalty::ThickPenalty)
    @unpack tmin, α, valid = thick_penalty
    if valid 
        tp0 = thickness_penalty(am0, tmin, α)
        tp1 = thickness_penalty(am1, tmin, α)
        return (tp1 - tp0)/δ
    else
        return 0.0
    end

    
end