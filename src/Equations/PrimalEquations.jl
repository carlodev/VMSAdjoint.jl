"""
Interface with SegregatedVMSSolver
"""



function equations_primal( simcase::Airfoil,params::Dict{Symbol,Any},time_dep::Symbol)
    equations_primal( simcase,params,Val(time_dep))
end

function equations_primal( simcase::Airfoil,params::Dict{Symbol,Any},::Val{:steady})
    @unpack uh = params

    _,_,Auu,Aup,Apu,App,_,_,rhs_v = segregated_equations(uh,params,simcase)
    res_prim((u, p), (v, q)) = Auu(u,v) + Aup(p,v) + Apu(u,q) + App(p,q)

    rhs((v,q)) = rhs_v(v)

    return res_prim, rhs
end



function equations_primal( simcase::Airfoil,params::Dict{Symbol,Any},::Val{:unsteady})
    @unpack dΩ = params
    @sunpack D = simcase
    Tuu,Tpu,Auu,Aup,Apu,App,_,_,rhs_v = segregated_equations(params[:uh],params,simcase)
    
    m(t, (u, p), (v, q)) = Tuu(u,v) + Tpu(u,q)

    res_prim(t,(u, p), (v, q)) = Auu(u,v) + Aup(p,v) + Apu(u,q) + App(p,q)

    rhs(t,(v,q)) = ∫(VectorValue(zeros(D)...) ⋅ v)dΩ

    return m, res_prim, rhs
end

