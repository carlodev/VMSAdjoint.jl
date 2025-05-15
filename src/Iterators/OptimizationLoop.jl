function solve_adjoint_optimization(opt)
   ## Set the optimization algorithm

    #Initialization
    w_init,_ = get_CST_values(des_points) 

    ###### Resolution
    (minf,minx,ret) = optimize(opt,w_init)
    
    return opt
end

function create_optimizer()
    ## Set the optimization algorithm
    opt = Opt(:LD_LBFGS, Ndes)
    ## Set the bounds
    opt.lower_bounds = [-0.1 .* ones(Int(Ndes/2));-1 .* ones(Int(Ndes/2))]
    opt.upper_bounds = [ones(Int(Ndes/2)); 0.1 .* ones(Int(Ndes/2))]

    ## Set the tolerance
    opt.xtol_rel = 2.5e-2

    opt.min_objective = optimization_loop

    return opt
end


"""
    optimization_loop(w::Vector,grad::Vector)

It receives:
-w: CST weights 
-grad: gradient dI/dw_i

It is using the solution at the previous iteration on the new geometry in order to reduce the initial transient
"""
function optimization_loop(w::Vector,grad::Vector)
 @unpack des_points,Ndes,shift,AoA, obj_fun, iter, 
            mesh_ref,uh00,ph00,UH,PH=params
    des_points = update_CST_weights(w,des_points)


    modelname =create_msh(des_points; AoA=AoA, mesh_ref=mesh_ref)
    model = GmshDiscreteModel(modelname)
    writevtk(model, "model_$iter")
    uh = uh00
    ph=ph00
  
    (uh,duhdt, UH, DUHDT), (ph, PH) = solve_inc_primal_u(model, params; filename="inc-results-$iter", uh00=uh00,ph00=ph00)
    
    
    
    #### Compute Unsteady CL
    CL_vec = zeros(length(PH))

    for (i,(uhval,phval)) in enumerate(zip(UH,PH))
    copyto!(uh.free_values, uhval)
    copyto!(ph.free_values, phval)
        fitnessval, CL  = IncompressibleAdjoint.obj_fun(model, params, uh,ph, compute_CL)
        CL_vec[i] = copy(CL)
    end

    params[:CL_history][iter+1] = copy(CL_vec)

    
    CD,CL = average_CD_CL(model,params,(uh,UH),(ph,PH))




    #### Adjoint Boundary Conditions
    params[:d_boundary] = VectorValue(0.0, (0.35-CL))
    

    ### run a steady adjoint on the average field
    average_field!(uh,UH[end-50:end])
    average_field!(ph,PH[end-50:end])

    #### average Cp
    airfoil_nodes_top, airfoil_nodes_bottom,_,_ = get_airfoil_characteristics(model, params; tag="airfoil")
    cp_top,cp_bottom = get_aerodynamic_features(params,model, uh,ph)
    params[:nodes_top][iter+1] = copy(getindex.(airfoil_nodes_top,1))
    params[:nodes_bottom][iter+1] =  copy(getindex.(airfoil_nodes_bottom,1))
    params[:cp_top][iter+1] =  copy(cp_top)
    params[:cp_bottom][iter+1] =  copy(cp_bottom)




    uhadj,phadj = solve_inc_adj_s(model, uh, ph,params; filename="res-adj-steady")

    J1ref,(J2_1ref,J2_2ref,J2_3ref,J2_4ref,J2_5ref)= IncompressibleAdjoint.compute_sensitivity(model,params, uh,ph,uhadj,phadj; objective_function=obj_fun)

    fitnessval, CL  = IncompressibleAdjoint.obj_fun(model, params, uh,ph, obj_fun)


    params[:iter]= iter +1 
    params[:βv][iter+1] = copy(w)

    params[:CL][iter+1] = copy(CL)
    params[:fitnessval][iter+1]= copy(fitnessval)
    J1=zeros(Ndes)
    J2_1=zeros(Ndes)
    J2_2=zeros(Ndes)
    J2_3=zeros(Ndes)
    J2_4=zeros(Ndes)
    J2_5=zeros(Ndes)

    copyto!(params[:uh00].free_values, uh.free_values)
    copyto!(params[:ph00].free_values, ph.free_values)


    for (i,ss) in enumerate(shift)

        des_temp = IncompressibleAdjoint.perturb_DesignParameters(des_points, i, ss)
        modelname_tmp =create_msh(des_temp; AoA=AoA, mesh_ref=mesh_ref)
        model_tmp = GmshDiscreteModel(modelname_tmp)

        J1tmp,(J2_1tmp,J2_2tmp,J2_3tmp,J2_4tmp,J2_5tmp)= IncompressibleAdjoint.compute_sensitivity(model_tmp,params, uh,ph,uhadj,phadj; objective_function=compute_CL)


        J1[i] = J1tmp
        J2_1[i] = J2_1tmp
        J2_2[i] = J2_2tmp
        J2_3[i] = J2_3tmp
        J2_4[i] = J2_4tmp
        J2_5[i] = J2_5tmp
        println(i)
    end
    
    J1s = (J1 .- J1ref)./δ
    J2s1 = (J2_1 .- J2_1ref)./δ
    J2s2 = (J2_2 .- J2_2ref)./δ
    J2s3 = (J2_3 .- J2_3ref)./δ
    J2s4 = (J2_4 .- J2_4ref)./δ
    J2s5 = (J2_5 .- J2_5ref)./δ
    
    
    J2s = J2s1 + J2s2 + J2s3 + J2s4 + J2s5
    Jtot = J1s+ J2s
    

    grad .= J2s.*des_bool

    params[:fgrad][iter+1] = copy(grad)
    println("CL = $CL")
    jldsave("results/NACA0012_1000_2p5/parameters_2.jld2"; params)
    return J1ref

end




