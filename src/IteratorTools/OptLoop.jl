

function solve_adjoint_optimization(adjp::AdjointProblem)
    @unpack solver = adjp

    #create initial
    w_init = get_DesignParameters(adjp.adesign)

    #create bounds
    Ndes = length(w_init)
    lb,ub = bounds_w(adjp.adesign, Ndes, solver.bounds)

    @info "Number of Design parameters: $Ndes"
    f, ∇f! = make_f_and_∇f(adjp, Ndes)


    opt_options = Optim.Options(iterations = solver.max_iter)  # change  to your desired limit

    ls = LineSearches.Static()
    step_options = LineSearches.InitialStatic(alpha=solver.αg, scaled=solver.scaled)
    # L-BFGS optimizer with line search control
    result = optimize(f, ∇f!,lb,ub, w_init, Fminbox(LBFGS(alphaguess=step_options, linesearch = ls)),opt_options)
    
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
    Nhalf = Int(Ndes ÷ 2)
    nose = 0.08
    lb = [nose; nose; nose; -0.15.* ones(Nhalf-3); -1.5 .* ones(Nhalf)]
    ub = [1.5 .* ones(Nhalf); nose; -nose;-nose; 0.15 .* ones(Nhalf-3)]
    return lb,ub
end


mutable struct SharedCache
    fval::Float64
    grad::Vector{Float64}
    valid::Bool
    iter::Int64
    uh
    ph
    adj_bc::Vector{Float64}
    am
    Cp
    adjp::AdjointProblem
    CDCL::Vector{Float64}
end

function make_f_and_∇f(adjp::AdjointProblem, N::Int64)
    cache = SharedCache(NaN, zeros(N), true, -1, nothing, nothing, [0.0,0.0], nothing,nothing,adjp, [0.0,0.0])

    x_last = similar(zeros(N))

    f(x) = begin
        if !isequal(x, x_last)
            copy!(x_last, x)
            cache.fval, cache.iter, cache.uh, cache.ph, cache.adj_bc, cache.am, cache.Cp, cache.CDCL = eval_f(x, cache)  # compute both
        end
        return cache.fval
    end

    ∇f!(g, x) = begin

        if !isequal(x, x_last)
            copy!(x_last, x)
            cache.fval, cache.iter, cache.uh, cache.ph, cache.adj_bc, cache.am, cache.Cp, cache.CDCL = eval_f(x, cache)
        end
        eval_∇f!(g, x, cache)  # must recompute if f not called
        copyto!(cache.grad, g)
    end

    return f, ∇f!
end



function eval_f(w::Vector, cache::SharedCache)

    @unpack  iter, uh, ph, adjp, am= cache
    @unpack JJfact, adesign, vbcase, timesol,solver = adjp

    meshinfo = vbcase.meshp.meshinfo
    physicalp = vbcase.simulationp.physicalp

    iter = iter+1
    @info "Iteration $iter started"
    println("Design parameters = $w")

    #create the new airfoil model from the weights w
    adesign = create_AirfoilDesign(adesign,w)
    modelname =create_msh(meshinfo,adesign, physicalp ; iter = iter)
    model = GmshDiscreteModel(modelname)

    writevtk(model, "model_$iter")
    am =  AirfoilModel(model, vbcase; am=am)


    #Solve Primal 
    if timesol==:steady
        filename = joinpath("Results_primal", "SOL_$(iter).vtu")
    else
           filename = "sol_$(iter)"
    end

    uh,ph = solve_inc_primal(am, vbcase, filename, timesol; uh0=uh,ph0=ph)    
 
    #### extract results from primal solution: Cp (and Cf)
    PressureCoefficient = get_aerodynamic_features(am,uh,ph)

    J,_ = JJfact

    
    fval, CDCL = obj_fun(am, vbcase, uh,ph, J)
    #### Adjoint Boundary Conditions
    #use the CLCD value to set the Boundary Condition
    adj_bc = -dJobj_fun(J, CDCL)
    

    return fval, iter, uh, ph, adj_bc, am, PressureCoefficient, CDCL
end


function eval_∇f!(grad::Vector, w::Vector,  cache::SharedCache)
    @unpack  iter, uh, ph, adjp, am, adj_bc, CDCL= cache
    @unpack JJfact,adesign, vbcase, timesol,solver = adjp
    
    _,Jfact = JJfact
    
    airfoil_case= vbcase
    adesign = create_AirfoilDesign(adesign,w)


    uhadj,phadj = solve_inc_adj(am, airfoil_case, adj_bc, "inc-adj-steady-$iter", :steady, uh, ph)

    Ndes = length(w) #number of design parameters
    δ = solver.δ #0.0001
    shift = CSTweights(Int(Ndes/2), δ)
    shiftv =   vcat(shift) #[δ,δ,δ,δ,δ...., -δ,-δ,-δ,-δ,.....]

    Jcorr = Jfact(CDCL) ## correction factor

    Jtot = iterate_perturbation(shiftv,adesign,am, airfoil_case,uh,uhadj, Jcorr )
    @info "Gradient: $Jtot"
    

    grad[:] = Jtot #update the gradients

    #update values iteration
    adj_sol = AdjSolution(iter, cache.fval, CDCL, w, grad, am, adj_bc, cache.Cp, uh, ph, uhadj, phadj  )

    mkpath("results")
    jldsave(joinpath("results","ADJ_SOL$(iter).jld2"); adj_sol)
    @info "Iteration $iter completed"


end


function iterate_perturbation(shift::Vector{Float64}, adesign::AirfoilDesign, am::AirfoilModel, airfoil_case::Airfoil, uh,uhadj, Jcorr::Float64 )
    Ndes = length(shift)

    meshinfo = airfoil_case.meshp.meshinfo
    physicalp =airfoil_case.simulationp.physicalp
    Ji = zeros(Ndes)
    for (i,ss) in enumerate(shift)
        @info "Perturbation Domain $i"

        adesign_tmp = perturb_DesignParameter(adesign, i, ss)
        modelname_tmp =create_msh(meshinfo,adesign_tmp, physicalp,"MeshPerturb"; iter = i+100)
        model_tmp = GmshDiscreteModel(modelname_tmp)
        am_tmp =  AirfoilModel(model_tmp, airfoil_case; am=am)

        Ji[i] = compute_sensitivity(am, am_tmp,ss, airfoil_case, uh,uhadj,Jcorr) 

    end
    
    return Ji
end
