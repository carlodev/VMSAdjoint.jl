using Revise
using AirfoilTools
using VMSAdjoint
using Gridap, GridapGmsh
using SegregatedVMSSolver.ParametersDef
using PartitionedArrays


fname = "n0012.csv" #airfoil coordinates to load
AoA = 4.0

ap0 = get_airfoil_coordinates(joinpath(@__DIR__, fname))

control_px = collect(LinRange(0.05,0.95,20))
control_points = ControlPoints(control_px,control_px)


R = 0.15 #support radius. If is bigger, the deformation radius increases, try increasing it
RBFfun = RBFFunctionLocalSupport(RBF_CP4, R)

rbfg = RBFGeometry(control_points,RBFfun)
rbfd = RBFDesign(rbfg, ap0)


sprob = StabilizedProblem(VMS(2))

physicalp = PhysicalParameters(Re=1000, u_in=[1.0,0.0])
timep = TimeParameters(dt=0.05, tF=15.0, time_window=(10.0, 15.0))

meshinfo = AirfoilMesh(AoA= AoA, meshref=2)
meshp = MeshParameters((2,2), 2, meshinfo)
exportp = ExportParameters(printinitial=true,printmodel=true,name_tags=["airfoil"], fieldexport=[["uh","ph","friction"]])

solverp = SolverParameters(θ=1.0)

simparams = SimulationParameters(timep,physicalp,solverp,exportp)

airfoil_case = Airfoil(meshp,simparams,sprob)


#Define the Objective Function, first argument [CD,CL], then you can define whatever you want

function J(CDCL; CLtarget=0.75)
    CD,CL=CDCL
    return -CL #0.5 * (CL - CLtarget)^2
end  



adj_solver = AdjSolver(δ=0.0001, opt_alg=:LD_LBFGS)

adjoint_airfoil_problem = AdjointProblem( rbfd,airfoil_case,adj_solver,:unsteady, J)



VMSAdjoint.solve(adjoint_airfoil_problem, with_debug)





# function get_nodes_idx(model::GridapDistributed.DistributedDiscreteModel, parts, D::Int64, AoA::Float64, tag::String; am=nothing)
#     Γ = BoundaryTriangulation(model; tags=tag)
#     n_Γ = get_normal_vector(Γ)

#     @assert D == 2 "only 2D supported"

#     if isnothing(params)
#         local_unique_idx = map(Γ.trians) do ttrian
#             visgrid = get_visgrid(ttrian)
#             visgrid_ = conv_VectorValue.(visgrid.sub_grid.node_coordinates)
#             unique_idx = unique(i -> visgrid_[i], eachindex(visgrid_)) #Indexes of unique nodes on each part
#             return unique_idx
#         end

#     else 

#         @unpack local_unique_idx = params

#     end

#     local_unique_nodes = map(Γ.trians,local_unique_idx) do ttrian, lui
#         visgrid = get_visgrid(ttrian)
#         visgrid_ = conv_VectorValue.(visgrid.sub_grid.node_coordinates)
#         nodes_tri = unwrap_vector(visgrid_[lui])
#         return nodes_tri
#     end

#     glun = gather(local_unique_nodes)

#     global_unique_idx = []
#     nodes_uniques = []
#     D = 2 #2D supported

#     if isnothing(params)
#     map(glun, parts) do g, part
#             if part == 1 #On Main procs
#                 gg = wrap_vector(g.data, D)
#                 global_unique_idx = unique(i -> gg[i], eachindex(gg))
#                 nodes_uniques = gg[global_unique_idx]
#             end
#         end
#     else
#         @unpack global_unique_idx = params
#         map(glun, parts) do g, part
#             if part == 1 #On Main procs
#                 gg = wrap_vector(g.data, D)
#                 nodes_uniques = gg[global_unique_idx]
#             end
#         end
  
#     end

#     airfoil_points = nodes_uniques
#     idx_uniques = global_unique_idx

#     if isnothing(params)
#         idx_above = is_above.(airfoil_points;AoA)
#         idx_top = findall(idx_above.> 0)
#         idx_bottom = findall(idx_above .< 0)
            
#         perm_top = sortperm(getindex.(airfoil_points[idx_top],1))
#         perm_bottom = sortperm(getindex.(airfoil_points[idx_bottom],1))
#         IDX_TOP = idx_top[perm_top]
#         IDX_BOTTOM = idx_bottom[perm_bottom]
#         IDX_TOP_UNIQUE=idx_uniques[IDX_TOP]
#         IDX_BOTTOM_UNIQUE=idx_uniques[IDX_BOTTOM]

#         params =Dict{Symbol,Any}(:local_unique_idx=> local_unique_idx, :global_unique_idx=>global_unique_idx)
#         params=merge!(params,Dict{Symbol,Any}(:IDX_TOP_UNIQUE=>IDX_TOP_UNIQUE,:IDX_BOTTOM_UNIQUE=>IDX_BOTTOM_UNIQUE,
#         :IDX_TOP=>IDX_TOP,:IDX_BOTTOM=>IDX_BOTTOM))

#         # merge!(params, Dict(:local_unique_idx=> local_unique_idx))
#         merge!(params, Dict(:nodes_uniques=> nodes_uniques))

#     else
        
#         @unpack IDX_TOP_UNIQUE,IDX_BOTTOM_UNIQUE,IDX_TOP,IDX_BOTTOM = params
#     end

#     # params=Dict(:IDX_TOP_UNIQUE=>IDX_TOP_UNIQUE,:IDX_BOTTOM_UNIQUE=>IDX_BOTTOM_UNIQUE,
#     # :IDX_TOP=>IDX_TOP,:IDX_BOTTOM=>IDX_BOTTOM)

#     return airfoil_points[IDX_TOP],airfoil_points[IDX_BOTTOM], params

# end


# function main_p(nparts,distribute)

#     parts  = distribute(LinearIndices((nparts,)))
#     modelname =create_msh(meshinfo,rbfd, physicalp)

#     model = GmshDiscreteModel(parts,modelname)
    
#     # writevtk(model, "model_p")
#     params = Dict(:parts=>parts)

#     AoA= 4.0
#     D = 2
#     nodesu,nodesl,params = get_nodes_idx(model,parts,D, AoA, "airfoil"; params=nothing)
#     nodesu,nodesl,params = get_nodes_idx(model,parts,D, AoA, "airfoil"; params=params)

#     ap = AirfoilPoints(getindex.(nodesu,1),getindex.(nodesl,1),getindex.(nodesu,2),getindex.(nodesl,2))
#     an = AirfoilNormals([zeros(2) for i =1:10] ,[zeros(2) for i =1:10] )
#     am = AirfoilModel(ap,model,an,params)
#     @unpack IDX_TOP_UNIQUE = am.params

#     reffeᵤ = ReferenceFE(lagrangian,VectorValue{2,Float64},1)
#     V = TestFESpace(model,reffeᵤ,conformity=:H1)
#     reffeₚ = ReferenceFE(lagrangian,Float64,1)
#     Q = TestFESpace(model,reffeₚ,conformity=:H1)

#     uh = interpolate_everywhere(VectorValue(rand(2)...), V)
#     ph = interpolate_everywhere(rand(), Q)

#     get_aerodynamic_featuresp(am, uh,ph, parts;tag="airfoil")

#     return am
# end


# am = with_debug() do distribute
#     main_p(nparts,distribute)
# end;






# function get_aerodynamic_featuresp(am, uh,ph;tag="airfoil")
#     @unpack IDX_TOP_UNIQUE,IDX_BOTTOM_UNIQUE, 
#     local_unique_idx,global_unique_idx,
#     IDX_TOP, IDX_BOTTOM = am.params
#     parts = am.parts

#     u_in = 1.0
#     q = 0.5 .* u_in^2

#     Γ = BoundaryTriangulation(am.model; tags=tag)

#     cf = Dict("uh"=>uh,"ph"=>ph)
   
#     fdat = GridapDistributed._prepare_fdata(Γ.trians, cf)

#     # pressure_top = pdata["ph"][IDX_TOP_UNIQUE]
#     # pressure_bottom  = pdata["ph"][IDX_BOTTOM_UNIQUE]

#     fieldh = map(Γ.trians, fdat,local_unique_idx) do ttrian, cf, unique_idx
#             visgrid = get_visgrid(ttrian)
#             pdata = Gridap.Visualization._prepare_pdata(ttrian, cf, visgrid.cell_to_refpoints)
#             field_h = pdata["ph"][unique_idx]
#             return field_h
#     end #end map

#     glun = gather(fieldh)
#     pressure_values_uniques=[]
#     map(glun, parts) do g, part
#         if part == 1 #On Main procs
#             # gg = wrap_vector(g.data, 2)
#             pressure_values_uniques = g.data[global_unique_idx]
#         end
#     end
#     println(pressure_values_uniques)
#     cp_top = pressure_values_uniques[IDX_TOP]./q
#     cp_bottom = pressure_values_uniques[IDX_BOTTOM]./q

#      return AirfoilScalar(cp_top,cp_bottom)
# end