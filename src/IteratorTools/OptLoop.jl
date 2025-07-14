

function solve_adjoint_optimization(adjp::AdjointProblem)
    @unpack solver = adjp
    @unpack ls, step_options = solver

    #create initial
    w_init = get_DesignParameters(adjp.adesign)

    #create bounds
    Ndes = length(w_init)
    lb,ub = bounds_w(adjp.adesign, Ndes, solver.bounds)

    @info "Number of Design parameters: $Ndes"
    f, ∇f! = make_f_and_∇f(adjp, Ndes)


    opt_options = Optim.Options(iterations = solver.max_iter)  # change  to your desired limit


    # L-BFGS optimizer with line search control
    #ls and step-options defined in the solver
    result = optimize(f, ∇f!,lb,ub, w_init, Fminbox(LBFGS(alphaguess=step_options, linesearch=ls)),opt_options)
    
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
    adesign
end

function make_f_and_∇f(adjp::AdjointProblem, N::Int64)
    cache = SharedCache(NaN, zeros(N), true, -1, nothing, nothing, [0.0,0.0], nothing,nothing,adjp, [0.0,0.0], nothing)

    x_last = similar(zeros(N))

    f(x) = begin
        if !isequal(x, x_last)
            copy!(x_last, x)
            cache.fval, cache.iter, cache.uh, cache.ph, cache.adj_bc, cache.am, cache.Cp, cache.CDCL, cache.adesign = eval_f(x, cache)  # compute both
        end
        return cache.fval
    end

    ∇f!(g, x) = begin

        if !isequal(x, x_last)
            copy!(x_last, x)
            cache.fval, cache.iter, cache.uh, cache.ph, cache.adj_bc, cache.am, cache.Cp, cache.CDCL,cache.adesign = eval_f(x, cache)
        end
        eval_∇f!(g, x, cache)  # must recompute if f not called
        copyto!(cache.grad, g)
    end

    return f, ∇f!
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


function eval_f(w::Vector, cache::SharedCache)

    @unpack  iter, uh, ph, adjp, am= cache
    @unpack J, adesign, vbcase, timesol,solver = adjp
    @unpack thick_penalty,regularization = solver

    meshinfo = vbcase.meshp.meshinfo
    physicalp = vbcase.simulationp.physicalp

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
    

    return fval, iter, uh, ph, adj_bc, am, PressureCoefficient, CDCL, adesign
end


function eval_∇f!(grad::Vector, w::Vector,  cache::SharedCache)
    @unpack  iter, uh, ph, adjp, am, adj_bc, CDCL, adesign= cache #take the updated adesign
    @unpack J, vbcase, timesol,solver = adjp
    @unpack thick_penalty = solver
    
    airfoil_case= vbcase
    
    uhadj,phadj = solve_inc_adj(am, airfoil_case, adj_bc, "inc-adj-"*string( timesol[2]) *"-$iter", timesol[2], uh, ph)

    Ndes = length(w) #number of design parameters
    δ = solver.δ #0.0001
    shift = CSTweights(Int(Ndes/2), δ)
    shiftv =   vcat(shift) #[δ,δ,δ,δ,δ...., -δ,-δ,-δ,-δ,.....]

    Ju,Jt = iterate_perturbation(shiftv,adesign,am, airfoil_case,solver, uh,uhadj )
    @info "Adjoint Gradient: $Ju"
    @info "Thickness Penalty Gradient: $Jt"


    grad[:] = Ju+Jt #update the gradients

    #update values iteration
    adj_sol = AdjSolution(iter, cache.fval, CDCL, w, grad, am, adj_bc, cache.Cp, uh, ph, uhadj, phadj  )

    mkpath("results")
    jldsave(joinpath("results","ADJ_SOL$(iter).jld2"); adj_sol)
    @info "Iteration $iter completed"


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

