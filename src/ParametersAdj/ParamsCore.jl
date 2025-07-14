
"""
    AirfoilNormals(nu,nl)

upper and lower sides normals
"""
struct AirfoilNormals
    nu::Vector{Vector{Float64}}
    nl::Vector{Vector{Float64}}
end

"""
    AirfoilModel

ap:: useful to plot them, and used to create the model
model:: is a Gridap object
params:: contains info on the nodes of the model, distinguishing from points on the top and bottom side
"""
struct AirfoilModel
    ap::AirfoilPoints
    model::Gridap.Geometry.UnstructuredDiscreteModel
    an::AirfoilNormals
    params::Dict{Symbol,Any}
end

"""
    AirfoilScalar

General structure for a scalar value on the surface of the airfoil (Cp or Cf)
"""
struct AirfoilScalar
    uval::Vector{Float64}
    lval::Vector{Float64}
end


function AirfoilScalar(am::AirfoilModel)
    nu = length(am.ap.xu)
    nl = length(am.ap.xl)
    return AirfoilScalar(zeros(nu), zeros(nl))
end

"""
    DesignBounds

Contains info on the boundary conditions for the design variables;
Check bounds_w to see how they are used
"""
@with_kw struct DesignBounds
    upper::Float64=0.5
    lower::Float64=-0.5
    Δy::Float64=0.005
end

"""
IdF(y) 
    Identity function - regularization not happening
"""
function IdF(x::Vector{Float64}, y::Vector{Float64}) 
    return y
end

@with_kw struct Regularization
    active::Bool = false
    iter_reg::Int64 = 2
    fun::Function = IdF
end

function thickness_penalty_default(pp; tmin=0.005, α=1_000.0)
    @assert α>0.0 "α has to be positive, α = $α not valid"

    xx, yu,yl = pp

    # Calculate thickness
    Δy = yu - yl

    # Check for minimum thickness violations
    thickness_violations = tmin .- Δy
    violations_locations = findall(thickness_violations.>0)
    !isempty(violations_locations) && println("Violations Locations x= $(xx[violations_locations])")

    #Σviolations^2
    penalty = sum((violation > 0) ? violation^2 : 0.0 for violation in thickness_violations)


    if isnan(penalty) || isinf(penalty)
        @error "Invalid penalty value detected. Check airfoil coordinates."
        return α   # Return a large penalty for invalid configurations
    end

    return penalty * α
end

"""
    ThickPenalty

Contains info on the thickness penalty constraint
"""
@with_kw struct ThickPenalty
    valid::Bool=true #decide to compute or not the penalty based on airfoil points
    thickness_penalty::Function= thickness_penalty_default
end


"""
    AdjSolver

Settings for the adjoint solver.

    Configurable options for opt_options
 x_abstol::Real = 0.0,
  x_reltol::Real = 0.0,
  f_abstol::Real = 0.0,
  f_reltol::Real = 0.0,
  g_abstol::Real = 1e-8,
  outer_x_abstol::Real = 0.0,
  outer_x_reltol::Real = 0.0,
  outer_f_abstol::Real = 0.0,
  outer_f_reltol::Real = 0.0,
  outer_g_abstol::Real = 1e-8,
  f_calls_limit::Int = 0,
  g_calls_limit::Int = 0,
  h_calls_limit::Int = 0,
  allow_f_increases::Bool = true,
  allow_outer_f_increases::Bool = true,
  successive_f_tol::Int = 1,
  iterations::Int = 1_000,
  outer_iterations::Int = 1000,
  store_trace::Bool = false,
  show_trace::Bool = false,
  extended_trace::Bool = false,
  show_warnings::Bool = true,
  show_every::Int = 1,
  callback = nothing,
  time_limit = NaN

"""
@with_kw struct AdjSolver
    δ::Float64=0.0001 #Perturbation of the design parameters for finite differences
    bounds::DesignBounds = DesignBounds()
    thick_penalty::ThickPenalty=ThickPenalty()
    regularization::Regularization=Regularization()
    ls = LineSearches.Static() #LineSearch
    step_options = LineSearches.InitialStatic(alpha=1.0, scaled=false)
    opt_options = Optim.Options()  # change  to your desired limit
end

@with_kw struct MeshSize
    BL_fl::Float64 = 1e-4 #fist layer height
    BL_tt::Float64 = 0.02 #total thickness
    airfoil_divisions::Int64 = 301
    meshref::Int64=1 #increase it, and it wil increase the resolution. Put = 0 and it will be very coarse, useful to debug
    H::Real = 8
    Lback::Real = 8
end

@with_kw struct AirfoilMesh <:MeshInfo
    AoA::Real #Angle of Attack - degrees
    folder::String="MeshFiles"
    MS::MeshSize = MeshSize()
end

struct AdjointProblem
    adesign::AirfoilDesign
    vbcase::VelocityBoundaryCase
    solver::AdjSolver
    timesol::Tuple{Symbol,Symbol}
    J::Function #(objective function)
    function AdjointProblem(   adesign::AirfoilDesign,vbcase::VelocityBoundaryCase,solver::AdjSolver,timesol::Tuple{Symbol,Symbol},J::Function)
        for at in timesol
            if at ∉ (:steady, :unsteady)
                throw(ArgumentError("Invalid timesol: $timesol. Must be :steady, :unsteady"))
            end
        end
        
        new(adesign, vbcase,solver,timesol,J)
    end
end

function AdjointProblem(    adesign::AirfoilDesign,vbcase::VelocityBoundaryCase,timesol::Tuple{Symbol,Symbol}, J::Function)
    solver=AdjSolver()
return     AdjointProblem(    adesign,vbcase,solver,timesol,J)
end

"""
    AdjSolution

It stores the results of the adjoint iteration.
fitnessval : value of the fitness function, it should go towards zero as the iterations are increasing
βv : value of the design parameters (e.g. CST weights) - it has been kept general
fgrad : ∂J/∂βv derivative of the objective function with respect to the design variable βv
airfoil_model : we store the arifoil model at this
"""
struct AdjSolution
    iter::Int64
    fitnessval::Float64 
    CDCL::Vector{Float64}
    βv::Vector{Float64} 
    fgrad::Vector{Float64} 
    airfoil_model::AirfoilModel
    bc_adj::Vector{Float64} 
    pressure_distribution::AirfoilScalar
    uh
    ph
    uhadj
    phadj
end




