

function solve_adjoint_optimization(adjp::AdjointProblem)
    ## Set the optimization algorithm
   
    opt,w_init = create_optimizer(adjp, optimization_loop)

    ###### Resolution
    (minf,minx,ret) = optimize(opt,w_init)
    
    return (minf,minx,ret) 
end

function create_optimizer(adjp::AdjointProblem, opt_loop::Function)
    @unpack solver = adjp

    #create initial
    w_init = get_DesignParameters(adjp.adesign)

    #create bounds
    Ndes = length(w_init)
    
    ## Set the optimization algorithm
    opt = Opt(solver.opt_alg, Ndes) #solver.opt_alg = :LD_LBFGS

    opt.initial_step = 0.01 .* ones(Ndes)  # or smaller

    lb,ub = bounds_w(adjp.adesign,  Ndes)

    opt.lower_bounds =lb
    opt.upper_bounds =ub
    

    ## Set the tolerance
    opt.xtol_rel = solver.tol

    #initialize loop parameters
    loop_params=  init_loop(adjp)

    opt_loop_obj = (w, grad) -> opt_loop(w, grad, loop_params)

    opt.min_objective = opt_loop_obj
    return opt,w_init
end



function bounds_w(adesign::RBFDesign, Ndes::Int64)
    Nhalf = Int(Ndes ÷ 2)
    lb = [-0.05 .* ones(Nhalf);-0.25.* ones(Nhalf)]
    ub =[0.25 .* ones(Nhalf); 0.05.* ones(Nhalf)]
    lb[1] = 0.005
    ub[Nhalf+1] = -0.005
    
    
    [@assert ub1>lb1 for (ub1,lb1) in zip(ub,lb)]

    return lb,ub
end


function bounds_w(adesign::AirfoilCSTDesign,  Ndes::Int64)
    Nhalf = Int(Ndes ÷ 2)
    nose = 0.08
    lb = [nose; nose; nose; -0.15.* ones(Nhalf-3); -1.5 .* ones(Nhalf)]
    ub = [1.5 .* ones(Nhalf); nose; -nose;-nose; 0.15 .* ones(Nhalf-3)]
    return lb,ub
end


function init_loop(adjp::AdjointProblem)
    loop_params = Dict{Symbol,Any}(:iter=> -1, :uh=>nothing, :ph=>nothing, :adjp=>adjp, :am=>nothing)
    return loop_params
end


"""
    optimization_loop(w::Vector,grad::Vector)

It receives:
-w: CST weights 
-grad: gradient dI/dw_i

It is using the solution at the previous iteration on the new geometry in order to reduce the initial transient
"""
function optimization_loop(w::Vector,grad::Vector, loop_params::Dict{Symbol,Any})

    @unpack  iter, uh, ph, adjp, am= loop_params
    @unpack J,adesign, vbcase, timesol,solver = adjp

    println("w = $w")

    #rename
    airfoil_case = vbcase


    meshinfo = airfoil_case.meshp.meshinfo
    physicalp = airfoil_case.simulationp.physicalp


    iter = iter+1

    iter>10 && throw(StopException("Maximum number of iterations"))


    @info "Iteration $iter started"

    #create the new airfoil model from the weights w
    adesign = create_AirfoilDesign(adesign,w)

    

    modelname =create_msh(meshinfo,adesign, physicalp ; iter = iter)
    model = GmshDiscreteModel(modelname)
    writevtk(model, "model_$iter")

    am =  AirfoilModel(model, airfoil_case; am=am)


    #Solve Primal 
    if timesol==:steady
        filename = joinpath("Results_primal", "SOL_$(iter).vtu")
    else
           filename = "sol_$(iter)"
    end

    uh,ph = solve_inc_primal(am, airfoil_case, filename, timesol; uh0=uh,ph0=ph)    
 
    #### extract results from primal solution: Cp (and Cf)
    PressureCoefficient = get_aerodynamic_features(am,uh,ph)


    fval, CDCL = obj_fun(am, airfoil_case, uh,ph, J)
 
    #### Adjoint Boundary Conditions
    #use the CLCD value to set the Boundary Condition
    adj_bc = -dJobj_fun(J, CDCL)

    uhadj,phadj = solve_inc_adj(am, airfoil_case, adj_bc, "inc-adj-steady-$iter", :steady, uh, ph)


    J1ref,(J2_1ref,J2_2ref,J2_3ref,J2_4ref,J2_5ref)= compute_sensitivity(airfoil_case, am,adj_bc, uh,ph,uhadj,phadj, J; i=iter)


    Ndes = length(w) #number of design parameters
    δ = solver.δ #0.01
    shift = CSTweights(Int(Ndes/2), δ)
    shiftv = vcat(shift)

    J1, (J2_1,J2_2,J2_3,J2_4,J2_5) = iterate_perturbation(shiftv,adesign,airfoil_case,adj_bc,uh,ph,uhadj,phadj, J  )
    J1s = (J1 .- J1ref)./shiftv
    J2s1 = (J2_1 .- J2_1ref)./shiftv
    J2s2 = (J2_2 .- J2_2ref)./shiftv
    J2s3 = (J2_3 .- J2_3ref)./shiftv
    J2s4 = (J2_4 .- J2_4ref)./shiftv
    J2s5 = (J2_5 .- J2_5ref)./shiftv
    
    
    J2s = J2s1 + J2s2 + J2s3 + J2s4 + J2s5
    Jtot = J1s .+ J2s
    

    grad .= Jtot


    #update values iteration


    adj_sol = AdjSolution(iter, fval, CDCL, w, grad, am, adj_bc, PressureCoefficient, uh,ph )

    loop_params[:iter] = iter
    loop_params[:uh] = uh
    loop_params[:ph] = ph
    loop_params[:am] = am

    jldsave("results/ADJ_SOL$(iter).jld2"; adj_sol)
    @info "Iteration $iter completed"

    return fval

end



function iterate_perturbation(shift::Vector{Float64},adesign::AirfoilDesign,airfoil_case::Airfoil,adj_bc::Vector{Float64}, uh,ph,uhadj,phadj, J::Function  )
    Ndes = length(shift)

    meshinfo = airfoil_case.meshp.meshinfo
    physicalp =airfoil_case.simulationp.physicalp
    
    J1=zeros(Ndes)    #Geometric term
    J2_1=zeros(Ndes)
    J2_2=zeros(Ndes)
    J2_3=zeros(Ndes)
    J2_4=zeros(Ndes)
    J2_5=zeros(Ndes)

    for (i,ss) in enumerate(shift)
        
        adesign_tmp = perturb_DesignParameter(adesign, i, ss)

        modelname_tmp =create_msh(meshinfo,adesign_tmp, physicalp,"MeshPerturb"; iter = i)
        model_tmp = GmshDiscreteModel(modelname_tmp)
        am_tmp =  AirfoilModel(model_tmp, airfoil_case)

        J1tmp,(J2_1tmp,J2_2tmp,J2_3tmp,J2_4tmp,J2_5tmp)= compute_sensitivity(airfoil_case, am_tmp,adj_bc, uh,ph,uhadj,phadj, J; i=i)

        J1[i] = J1tmp
        J2_1[i] = J2_1tmp
        J2_2[i] = J2_2tmp
        J2_3[i] = J2_3tmp
        J2_4[i] = J2_4tmp
        J2_5[i] = J2_5tmp
        println("Perturbation Domain $i")
    end
    
    return J1, (J2_1,J2_2,J2_3,J2_4,J2_5)
end
