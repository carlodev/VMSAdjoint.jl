

function solve_adjoint_optimization(adjp::AdjointProblem)
    @unpack solver = adjp
    @unpack opt_alg = solver

    #create initial
    w_init = get_DesignParameters(adjp.adesign)

    #create bounds
    Ndes = length(w_init)
    lb,ub = bounds_w(adjp.adesign, Ndes, solver.bounds)

    [@assert ub1 >= w1 >= lb1 "Parameter number:$(i) not valid: $(ub1) $(w1) $(lb1) not in-bound" for (i,(ub1,w1,lb1)) in enumerate(zip(ub,w_init,lb)) ]

    @info "Number of Design parameters: $Ndes"
    opt = NLopt.Opt(opt_alg, Ndes)
    NLopt.lower_bounds!(opt, lb)
    NLopt.upper_bounds!(opt, ub)
    
    scache = SharedCache(-1,nothing,nothing, adjp)

    opt_loop_obj = (x,g) -> f_and_∇f(x,g, scache)
    
    # L-BFGS optimizer with line search control
    #ls and step-options defined in the solver
    NLopt.min_objective!(opt, opt_loop_obj)

    min_f, min_x, ret = NLopt.optimize(opt, w_init)

    
    return true
end



#Boundary conditions on design parameters

function bounds_w(adesign::RBFDesign, Ndes::Int64,bounds::DesignBounds )
    @unpack upper, lower, Δy = bounds
    Nhalf = Int(Ndes ÷ 2)
    lb = [-Δy .* ones(Nhalf);lower.* ones(Nhalf)]
    ub =[upper .* ones(Nhalf); Δy.* ones(Nhalf)]
    lb[1] = Δy
    ub[Nhalf+1] = -Δy
    
    @info "Bounds settings"
    println("lower:$lb")
    println("upper:$ub")
 
    [@assert ub1>lb1 for (ub1,lb1) in zip(ub,lb)]

    return lb,ub
end


function bounds_w(adesign::AirfoilCSTDesign,  Ndes::Int64, bounds::DesignBounds)
    @unpack upper, lower, Δy = bounds

    Nhalf = Int(Ndes ÷ 2)
    nose = 0.08
    lb = [nose; nose; nose; -0.15.* ones(Nhalf-3); -1.5 .* ones(Nhalf)]
    ub = [1.5 .* ones(Nhalf); nose; -nose;-nose; 0.15 .* ones(Nhalf-3)]
    
    @info "Bounds settings"
    println("lower:$lb")
    println("upper:$ub")

    return lb,ub
end


mutable struct SharedCache
    iter::Int64
    uh
    ph
    adjp::AdjointProblem
end

function f_and_∇f(w::Vector{Float64}, grad::Vector{Float64}, cache::SharedCache)

    @unpack iter, uh, ph, adjp = cache
    @unpack J, adesign, vbcase, timesol,solver = adjp
    @unpack thick_penalty,regularization = solver

    meshinfo = vbcase.meshp.meshinfo
    physicalp = vbcase.simulationp.physicalp

    Ndes = length(w)
    
    iter = iter+1
    @info "Iteration $iter started"
    println("Design parameters = $w")


    #create the new airfoil model from the weights w
    adesign = create_AirfoilDesign(adesign,w)
    adesign = regularize_airfoil(adesign, iter, regularization) #design regularization
    model = generate_regularized_model(adesign, iter, 0.0, meshinfo, physicalp, "MeshFiles")
    
    writevtk(model, "model_$iter")
    am =  AirfoilModel(model, vbcase)


    #Solve Primal 
    if timesol==:steady
        filename = joinpath("Results_primal", "SOL_$(iter).vtu")
    else
        filename = "sol_$(iter)"
    end

    uh,ph = solve_inc_primal(am, vbcase, filename, timesol[1]; uh0=nothing,ph0=nothing)    
 
    #### extract results from primal solution: Cp (and Cf)
    PressureCoefficient = get_aerodynamic_features(am,uh,ph)
    
    fval, CDCL = obj_fun(am, vbcase, uh,ph, thick_penalty, J)
    #### Adjoint Boundary Conditions
    #use the CLCD value to set the Boundary Condition
    
    adj_bc = -dJobj_fun(J, CDCL)
    
    #### Start Computing the Gradient ####

    uhadj,phadj = solve_inc_adj(am, vbcase, adj_bc, "inc-adj-"*string( timesol[2]) *"-$iter", timesol[2], uh, ph)

    δ = solver.δ #0.0001
    shift = CSTweights(Int(Ndes/2), δ)
    shiftv =   vcat(shift) #[δ,δ,δ,δ,δ...., -δ,-δ,-δ,-δ,.....]

    Ju,Jt = iterate_perturbation(shiftv,adesign,am, vbcase,solver, uh,uhadj )
    @info "Adjoint Gradient: $Ju"
    @info "Thickness Penalty Gradient: $Jt"
    grad[:] = Ju+Jt #update the gradients


    ### Save the Solution ###
    #update values iteration
    adj_sol = AdjSolution(iter, fval, CDCL, w, grad, am, adj_bc, PressureCoefficient, uh, ph, uhadj, phadj  )

    mkpath("results")
    jldsave(joinpath("results","ADJ_SOL$(iter).jld2"); adj_sol)
    @info "Iteration $iter completed"

    ### update the cache
    cache.iter = iter
    cache.uh = uh
    cache.ph = ph



    return fval
end



function iterate_perturbation(shift::Vector{Float64}, adesign::AirfoilDesign, am::AirfoilModel, airfoil_case::Airfoil, solver::AdjSolver, uh,uhadj )
    Ndes = length(shift)
    @unpack thick_penalty,regularization = solver

    meshinfo = airfoil_case.meshp.meshinfo
    physicalp =airfoil_case.simulationp.physicalp
    Ji = zeros(Ndes)
    Jthickness = zeros(Ndes)
    for (i,ss) in enumerate(shift)
        @info "Perturbation Domain $i"

        model_tmp = generate_regularized_model(adesign, i, ss, meshinfo, physicalp, "MeshPerturb")
        am_tmp =  AirfoilModel(model_tmp, airfoil_case)

        Ji[i],Jthickness[i] = compute_sensitivity(am, am_tmp,adesign, i,ss, airfoil_case,thick_penalty, uh,uhadj) 

    end
    
    return Ji, Jthickness
end

