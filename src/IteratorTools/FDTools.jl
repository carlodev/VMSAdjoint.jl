"""
    finite_difference_analysis(adjoint_airfoil_problem::AdjointProblem)
Iterator for Finite Difference Analysis
"""
function finite_difference_analysis(adjoint_airfoil_problem::AdjointProblem; idxs::Vector{Int64}=Int64[])
    @unpack timesol, adesign, J, vbcase, solver= adjoint_airfoil_problem
    @unpack thick_penalty = solver
    @sunpack order = vbcase


    fddir = "FD"
    mkpath(fddir)

    airfoil_case = vbcase
    w = get_DesignParameters(adesign)

    meshinfo = airfoil_case.meshp.meshinfo
    physicalp =airfoil_case.simulationp.physicalp

    Ndes = length(w) #number of design parameters

    if isempty(idxs)
        idxs = collect(1:1:Ndes)
    else
        @assert maximum(idxs) <= Ndes "The IDX of Finite Difference has to be in the range of tha maximum number of design parameters; maximum(idxs) = $(maximum(idxs)), Number of design parameters: $Ndes"
    end

    
    δ = solver.δ #0.001
    shift = CSTweights(Int(Ndes/2), δ)
    shiftv =   vcat(shift)
    
    #create the new airfoil model from the weights w
    adesign = create_AirfoilDesign(adesign,w)
    modelname =create_msh(meshinfo,adesign, physicalp ; iter = 0)
    model = GmshDiscreteModel(modelname)
    writevtk(model, joinpath(fddir,"model_FD_0"))
    am =  AirfoilModel(model, airfoil_case; am=nothing)

    filename = joinpath(fddir, "FD_0")

   

    uh0,ph0 = solve_inc_primal(am, airfoil_case, filename, timesol)    
    fval0, CLCD0 = obj_fun(am, airfoil_case, uh0,ph0,thick_penalty, J)
    fval_fd, CLCD_fd = iterate_fd(shiftv,idxs, adesign,am,airfoil_case,timesol,thick_penalty, J )

    fval_grad = (fval_fd .- fval0)./shiftv
    CLCD_grad= map(CLCDi-> CLCDi .- CLCD0, CLCD_fd) ./shiftv
    
    sol_fd = Dict(:CLCD_grad=>CLCD_grad, :fval_grad=>fval_grad, :idxs=>idxs)
    jldsave( joinpath(fddir, "FD_SOL.jld2"); sol_fd)


    return fval_grad, CLCD_grad

end

function iterate_fd(shift::Vector{Float64}, idxs::Vector{Int64}, adesign::AirfoilDesign,am::AirfoilModel,airfoil_case::Airfoil,timesol::Symbol,thick_penalty, J::Function  )
    fddir = "FD"
    mkpath(fddir)

    # Ndes = length(shift)
    Nidxs = length(idxs)

    meshinfo = airfoil_case.meshp.meshinfo
    physicalp =airfoil_case.simulationp.physicalp


    fval_fd = zeros(Nidxs)
    CLCD_fd = [zeros(length(2)) for i= 1:Nidxs]

    
    for i in idxs
        
        ss = shift[i]

        println("Perturbation Domain $i")

        adesign_tmp = perturb_DesignParameter(adesign, i, ss)

        modelname_tmp =create_msh(meshinfo,adesign_tmp, physicalp,"MeshPerturb"; iter = i+100)
        model_tmp = GmshDiscreteModel(modelname_tmp)
        am_tmp =  AirfoilModel(model_tmp, airfoil_case, am=am)

        filename = joinpath(fddir, "FD_$(i)")

        
        uh_tmp,ph_tmp = solve_inc_primal(am_tmp, airfoil_case, filename, timesol; uh0=nothing,ph0=nothing)    

        fval_fd[i], CLCD_fd[i] = obj_fun(am_tmp, airfoil_case, uh_tmp,ph_tmp,thick_penalty, J)

    end
    
    return fval_fd, CLCD_fd
end

