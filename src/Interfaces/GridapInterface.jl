
# uniqueidx(v) = unique(i -> v[i], eachindex(v))

function uniqueidx(v::AbstractVector)
    tol = 1e-7
    idxs = Int[]
    for (i, vi) in enumerate(v)
        if all(j -> norm(vi - v[j]) > tol, idxs)
            push!(idxs, i)
        end
    end
    return idxs
end



function get_nodes_idx(model::GridapDistributed.DistributedDiscreteModel, parts, D::Int64, AoA::Float64, tag::String; params=nothing)
    Γ = BoundaryTriangulation(model; tags=tag)

    println("In get_nodes_idx")

    @assert D == 2 "only 2D supported"

    if isnothing(params)
        local_unique_idx = map(Γ.trians) do ttrian
            visgrid = get_visgrid(ttrian)
            visgrid_ = conv_VectorValue.(visgrid.sub_grid.node_coordinates)
            unique_idx = unique(i -> visgrid_[i], eachindex(visgrid_)) #Indexes of unique nodes on each part
            return unique_idx
        end

    else 

        @unpack local_unique_idx = params

    end

    local_unique_nodes = map(Γ.trians,local_unique_idx) do ttrian, lui
        visgrid = get_visgrid(ttrian)
        visgrid_ = conv_VectorValue.(visgrid.sub_grid.node_coordinates)
        nodes_tri = unwrap_vector(visgrid_[lui])
        return nodes_tri
    end

    glun = gather(local_unique_nodes)

    global_unique_idx = []
    nodes_uniques = []
    D = 2 #2D supported

    if isnothing(params)
    map(glun, parts) do g, part
            if part == 1 #On Main procs
                gg = wrap_vector(g.data, D)
                global_unique_idx = unique(i -> gg[i], eachindex(gg))
                nodes_uniques = gg[global_unique_idx]
            end
        end
    else
        @unpack global_unique_idx = params
        map(glun, parts) do g, part
            if part == 1 #On Main procs
                gg = wrap_vector(g.data, D)
                nodes_uniques = gg[global_unique_idx]
            end
        end
  
    end

    airfoil_points = nodes_uniques
    idx_uniques = global_unique_idx

    if isnothing(params)
        idx_above = is_above.(airfoil_points;AoA)
        idx_top = findall(idx_above.> 0)
        idx_bottom = findall(idx_above .< 0)
            
        perm_top = sortperm(getindex.(airfoil_points[idx_top],1))
        perm_bottom = sortperm(getindex.(airfoil_points[idx_bottom],1))
        IDX_TOP = idx_top[perm_top]
        IDX_BOTTOM = idx_bottom[perm_bottom]
        IDX_TOP_UNIQUE=idx_uniques[IDX_TOP]
        IDX_BOTTOM_UNIQUE=idx_uniques[IDX_BOTTOM]

        params =Dict{Symbol,Any}(:local_unique_idx=> local_unique_idx, :global_unique_idx=>global_unique_idx)
        params=merge!(params,Dict{Symbol,Any}(:IDX_TOP_UNIQUE=>IDX_TOP_UNIQUE,:IDX_BOTTOM_UNIQUE=>IDX_BOTTOM_UNIQUE,
        :IDX_TOP=>IDX_TOP,:IDX_BOTTOM=>IDX_BOTTOM))

        # merge!(params, Dict(:local_unique_idx=> local_unique_idx))
        merge!(params, Dict(:nodes_uniques=> nodes_uniques))

    else
        
        @unpack IDX_TOP_UNIQUE,IDX_BOTTOM_UNIQUE,IDX_TOP,IDX_BOTTOM = params
    end

    # params=Dict(:IDX_TOP_UNIQUE=>IDX_TOP_UNIQUE,:IDX_BOTTOM_UNIQUE=>IDX_BOTTOM_UNIQUE,
    # :IDX_TOP=>IDX_TOP,:IDX_BOTTOM=>IDX_BOTTOM)

    return airfoil_points[IDX_TOP],airfoil_points[IDX_BOTTOM], params

end



# function get_nodes_idx(model, AoA::Float64, tag::String; params=nothing)

#     f = (reffe) -> Gridap.Geometry.UnstructuredGrid(reffe)
#     Γ = BoundaryTriangulation(model; tags=tag)

#     trian_tag = Γ.trian
#     ref_grids0 = map(f, Gridap.Geometry.get_reffes(trian_tag))
#     visgrid0 = Gridap.Visualization.VisualizationGrid(trian_tag, ref_grids0)
#     airfoil_points0 = visgrid0.sub_grid.node_coordinates
    
#     idx_uniques = uniqueidx(airfoil_points0)
#     airfoil_points = airfoil_points0[idx_uniques]
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
#     else
#         @unpack IDX_TOP_UNIQUE,IDX_BOTTOM_UNIQUE,IDX_TOP,IDX_BOTTOM = params
#     end

#     params=Dict(:IDX_TOP_UNIQUE=>IDX_TOP_UNIQUE,:IDX_BOTTOM_UNIQUE=>IDX_BOTTOM_UNIQUE,
#     :IDX_TOP=>IDX_TOP,:IDX_BOTTOM=>IDX_BOTTOM)


#     return airfoil_points[IDX_TOP],airfoil_points[IDX_BOTTOM], params
# end




# function get_normals(model, params, tag)

#     @unpack IDX_TOP, IDX_BOTTOM = params
#     Γ = BoundaryTriangulation(model; tags=tag)
#     n_Γ = -1 .*get_normal_vector(Γ) #pointing outward
#     n_Γ_dict = Dict("n_Γ" => n_Γ)
#     f = (reffe) -> Gridap.Geometry.UnstructuredGrid(reffe)

#     n_trian_tag = n_Γ.trian
#     n_ref_grids0 = map(f, Gridap.Geometry.get_reffes(n_trian_tag))
#     n_visgrid0 = Gridap.Visualization.VisualizationGrid(n_trian_tag, n_ref_grids0)
#     n_Γ_data = Gridap.Visualization._prepare_pdata(n_trian_tag, n_Γ_dict, n_visgrid0.cell_to_refpoints)
#     airfoil_normals = n_Γ_data["n_Γ"]
    

#     airfoil_normals_top = airfoil_normals[IDX_TOP]
#     airfoil_normals_bottom = airfoil_normals[IDX_BOTTOM]
#     return airfoil_normals_top,airfoil_normals_bottom
# end

function ParametersAdj.AirfoilModel(model, mcase::Airfoil, parts, am::Nothing; tag="airfoil")

    println("In AM")
    @sunpack AoA,D = mcase
    nodesu,nodesl,params = get_nodes_idx(model, parts, D, AoA, tag;)



    ap = AirfoilPoints(getindex.(nodesu,1),getindex.(nodesl,1),getindex.(nodesu,2),getindex.(nodesl,2))
    # nu,nl = get_normals(model, params, tag)
    # an = AirfoilNormals(map(n-> [n...], nu),map(n-> [n...], nl))
    an = AirfoilNormals([zeros(2) for i =1:10] ,[zeros(2) for i =1:10] )

    return AirfoilModel(ap,model,an,params, parts)
end



function ParametersAdj.AirfoilModel(model, mcase::Airfoil, am::AirfoilModel; tag="airfoil")
    @sunpack AoA,D = mcase
    nodesu,nodesl,params =     get_nodes_idx(model, am.parts, D, AoA, tag; params=am.params)
    ap = AirfoilPoints(getindex.(nodesu,1),getindex.(nodesl,1),getindex.(nodesu,2),getindex.(nodesl,2))
    # nu,nl = get_normals(model, params, tag)
    # an = AirfoilNormals(map(n-> [n...], nu),map(n-> [n...], nl))
    an = AirfoilNormals([zeros(2) for i =1:10] ,[zeros(2) for i =1:10] )

    return AirfoilModel(ap,model,an,params,am.parts)
end



function get_aerodynamic_features(am, uh,ph;tag="airfoil")
    @unpack IDX_TOP_UNIQUE,IDX_BOTTOM_UNIQUE, 
    local_unique_idx,global_unique_idx,
    IDX_TOP, IDX_BOTTOM = am.params
    parts = am.parts
    
    u_in = 1.0
    q = 0.5 .* u_in^2

    Γ = BoundaryTriangulation(am.model; tags=tag)

    cf = Dict("uh"=>uh,"ph"=>ph)
   
    fdat = GridapDistributed._prepare_fdata(Γ.trians, cf)

    # pressure_top = pdata["ph"][IDX_TOP_UNIQUE]
    # pressure_bottom  = pdata["ph"][IDX_BOTTOM_UNIQUE]

    fieldh = map(Γ.trians, fdat,local_unique_idx) do ttrian, cf, unique_idx
            visgrid = get_visgrid(ttrian)
            pdata = Gridap.Visualization._prepare_pdata(ttrian, cf, visgrid.cell_to_refpoints)
            field_h = pdata["ph"][unique_idx]
            return field_h
    end #end map

    glun = gather(fieldh)
    pressure_values_uniques=[]
    map(glun, parts) do g, part
        if part == 1 #On Main procs
            # gg = wrap_vector(g.data, 2)
            pressure_values_uniques = g.data[global_unique_idx]
        end
    end
    println(pressure_values_uniques)
    cp_top = pressure_values_uniques[IDX_TOP]./q
    cp_bottom = pressure_values_uniques[IDX_BOTTOM]./q

     return AirfoilScalar(cp_top,cp_bottom)
end



# function get_aerodynamic_features(am::AirfoilModel, uh,ph;tag="airfoil")

#     @unpack IDX_TOP_UNIQUE,IDX_BOTTOM_UNIQUE=am.params
#     u_in = 1.0
#     q = 0.5 .* u_in^2

#     f = (reffe) -> Gridap.Geometry.UnstructuredGrid(reffe)

#     Γ = BoundaryTriangulation(am.model; tags=tag)
#     cf = Dict("uh"=>uh,"ph"=>ph)
#     trian_tag = Γ.trian
    
#     ref_grids0 = map(f, Gridap.Geometry.get_reffes(trian_tag))
#     visgrid0 = Gridap.Visualization.VisualizationGrid(trian_tag, ref_grids0)
    

    
#     pdata = Gridap.Visualization._prepare_pdata(Γ,cf,visgrid0.cell_to_refpoints)

    
#     pressure_top = pdata["ph"][IDX_TOP_UNIQUE]
#     pressure_bottom  = pdata["ph"][IDX_BOTTOM_UNIQUE]
#     cp_top = pressure_top ./ q
#     cp_bottom = pressure_bottom./ q

#      return AirfoilScalar(cp_top,cp_bottom)
# end


function is_above(p; AoA)

    p2 = Point(cosd(AoA),-sind(AoA))
    return is_above(p, p2)
end

function is_above(p, p2; der_slope=1.0)
    der = p2[2]/p2[1] .* der_slope # -0.20 #derivative in trailing edge Manage the (-0.20) factor
    c = p2[2]/p2[1] #-der/8 +  #derivative in leading edge
    
    a = (2*p2[2]-p2[1]*der - c*p2[1])/(-p2[1]^3)
    
    b =(der-c-3*a*p2[1]^2)/(2 * p2[1])
    
    treshold_fun(x) = a*x^3 +b*x^2 +c*x
    
        res = p[2] - treshold_fun(p[1])
        if res > 0
            return 1
        else
            return -1
        end
end

