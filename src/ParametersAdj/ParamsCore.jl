abstract type MeshStructure end
struct Structured   <: MeshStructure end
struct Unstructured <: MeshStructure end


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
    model::Gridap.Geometry.DiscreteModel
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

Bounds on the design variables; consumed by `bounds_w` (one method per
parametrization).

- RBF designs use `upper`, `lower`, `Δy`.
- CST designs use `upper`/`lower` as the main weight bounds plus:
  `nose` — minimum magnitude of the three leading-edge weights (keeps a
  round nose), and `cross` — how far the remaining weights of one side may
  cross zero towards the other side.
"""
@with_kw struct DesignBounds
    upper::Float64=0.5
    lower::Float64=-0.5
    Δy::Float64=0.005
    nose::Float64=0.08
    cross::Float64=0.15
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

Options of the adjoint optimization loop.

`linear_backend` selects how the linear systems of the primal/adjoint solvers
are solved:
- `:cpu`   UMFPACK LU (default)
- `:gpu`   NVIDIA cuDSS LU; requires `using CUDA, CUDSS` before solving
- `:petsc` PETSc KSP driven by `SolverParameters.petsc_options`; the whole run
           must be wrapped in `with_petsc(simcase) do ... end`
- `:block` GridapSolvers FGMRES with a block velocity/pressure preconditioner

`mesh_motion` selects how the perturbed geometries of the sensitivity loop are
obtained:
- `true`  (default) deform the baseline mesh — no Gmsh call per design
          parameter, guaranteed topology match, exact mesh velocity
- `false` regenerate the mesh with Gmsh for every perturbation (previous
          behaviour; the fallback when a design step is too large to deform)

`mesh_stiffening` is the exponent χ of the cell-wise mesh stiffness
`(V_max/V_cell)^χ` used by the deformation; raise it if boundary-layer cells
distort, lower it towards 0 for a plain harmonic extension.

`adj_transient_steps` is how many time steps of adjoint start-up transient the
unsteady adjoint discards from its averaging (see
[`IncompressibleSolvers.solve_inc_adj_unsteady`](@ref)). It is a convergence
parameter, not a physical constant: sweep it, together with the primal
averaging window, before trusting an unsteady gradient.
"""
@with_kw struct AdjSolver
    δ::Float64=0.0001 #Perturbation of the design parameters for finite differences
    bounds::DesignBounds = DesignBounds()
    thick_penalty::ThickPenalty=ThickPenalty()
    regularization::Regularization=Regularization()
    opt_alg::Symbol = :LD_LBFGS
    linear_backend::Symbol = :cpu
    mesh_motion::Bool = true
    mesh_stiffening::Float64 = 1.0
    adj_transient_steps::Int = 10
    @assert linear_backend in (:cpu, :gpu, :petsc, :block) "linear_backend must be one of (:cpu, :gpu, :petsc, :block)"
    @assert mesh_stiffening >= 0.0 "mesh_stiffening (χ) must be non-negative"
    @assert adj_transient_steps >= 0 "adj_transient_steps must be non-negative"
end

@with_kw struct MeshSize
    BL_fl::Float64 = 1e-4 #fist layer height
    BL_tt::Float64 = 0.02 #total thickness
    airfoil_divisions::Int64 = 301
    meshref::Real=1.0 #increase it, and it wil increase the resolution. Put = 0 and it will be very coarse, useful to debug
    H::Real = 8
    Lback::Real = 8
    curv_turn_deg::Float64 = 5.0 #max turning angle (deg) of one airfoil surface element; 0 ⇒ uniform sizing
    @assert curv_turn_deg >= 0.0 "curv_turn_deg must be non-negative (0 disables curvature-based sizing)"
end

"""
    AirfoilMesh{Structured|Unstructured}

Mesh recipe for the airfoil domain. Supported element types:
- `Structured`:   `:TRI` or `:QUAD` (transfinite C-type mesh, optionally recombined)
- `Unstructured`: `:TRI` only (a quad/mixed unstructured mesh cannot be read by GridapGmsh)

The combination is validated in `create_msh`; CST designs additionally require a
`Structured` mesh (see `validate_mesh_config`).
"""
@with_kw struct AirfoilMesh{S<:MeshStructure} <:MeshInfo
    elements::Symbol=:TRI
    AoA::Real #Angle of Attack - degrees
    folder::String="MeshFiles"
    MS::MeshSize = MeshSize()
end

# AirfoilMesh(; kwargs...) = AirfoilMesh{Structured}(; kwargs...) #Default Constructor for AirfoilMesh


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




