using AirfoilTools
using Gridap
using Parameters
using SegregatedVMSSolver



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



function remove_outliers_from_curve(points::Vector, tol::Float64)

    filtered_pts = Vector{typeof(points[1])}()
    idxs = [1]
    push!(filtered_pts, points[1])  # Keep first point

    for i in 2:length(points)-1
        p_prev = points[i-1]
        p = points[i]
        p_next = points[i+1]
        
        # Expected position from linear interpolation between neighbors
        expected = 0.5 .* (p_prev + p_next)
        deviation = norm(p - expected)        
        if deviation < tol
            push!(filtered_pts, p)
            push!(idxs, i)

        end
    end

    push!(filtered_pts, points[end])  # Keep last point
    push!(idxs, length(points))

    return filtered_pts, idxs
end


function get_nodes_idx(model, AoA::Float64, tag::String)

    f = (reffe) -> Gridap.Geometry.UnstructuredGrid(reffe)
    Γ = BoundaryTriangulation(model; tags=tag)
    trian_tag = Γ.trian
    ref_grids0 = map(f, Gridap.Geometry.get_reffes(trian_tag))
    visgrid0 = Gridap.Visualization.VisualizationGrid(trian_tag, ref_grids0)
    airfoil_points0 = visgrid0.sub_grid.node_coordinates
    idx_uniques = uniqueidx(airfoil_points0)


    n_Γ = -1 .*get_normal_vector(Γ) #pointing outward
    n_Γ_dict = Dict("n_Γ" => n_Γ)
    f = (reffe) -> Gridap.Geometry.UnstructuredGrid(reffe)
    
    n_trian_tag = n_Γ.trian
    n_ref_grids0 = map(f, Gridap.Geometry.get_reffes(n_trian_tag))
    n_visgrid0 = Gridap.Visualization.VisualizationGrid(n_trian_tag, n_ref_grids0)
    n_Γ_data = Gridap.Visualization._prepare_pdata(n_trian_tag, n_Γ_dict, n_visgrid0.cell_to_refpoints)["n_Γ"]

    
    airfoil_points = airfoil_points0[idx_uniques]
    airfoil_normals =n_Γ_data[idx_uniques]
        
    idx_top = findall( map(a-> a[2] > 0.1 , airfoil_normals))
    idx_bottom = findall( map(a-> a[2] < -0.1 , airfoil_normals))

    perm_top = sortperm(getindex.(airfoil_points[idx_top],1))
    perm_bottom = sortperm(getindex.(airfoil_points[idx_bottom],1))
    IDX_TOP = idx_top[perm_top]
    IDX_BOTTOM = idx_bottom[perm_bottom]
    IDX_TOP_UNIQUE=idx_uniques[IDX_TOP]
    IDX_BOTTOM_UNIQUE=idx_uniques[IDX_BOTTOM]

    airfoil_points[IDX_TOP],airfoil_points[IDX_BOTTOM]

    airfoil_points_top,idx_top = remove_outliers_from_curve(airfoil_points[IDX_TOP], 0.01)
    airfoil_points_bottom,idx_bottom = remove_outliers_from_curve(airfoil_points[IDX_BOTTOM], 0.01)

    params=Dict(:IDX_TOP_UNIQUE=>IDX_TOP_UNIQUE[idx_top],:IDX_BOTTOM_UNIQUE=>IDX_BOTTOM_UNIQUE[idx_bottom],
    :IDX_TOP=>IDX_TOP,:IDX_BOTTOM=>IDX_BOTTOM)


    return airfoil_points_top,airfoil_points_bottom, params
end




function get_normals(model, params, tag)

    @unpack IDX_TOP, IDX_BOTTOM = params
    Γ = BoundaryTriangulation(model; tags=tag)
    n_Γ = -1 .*get_normal_vector(Γ) #pointing outward
    n_Γ_dict = Dict("n_Γ" => n_Γ)
    f = (reffe) -> Gridap.Geometry.UnstructuredGrid(reffe)

    n_trian_tag = n_Γ.trian
    n_ref_grids0 = map(f, Gridap.Geometry.get_reffes(n_trian_tag))
    n_visgrid0 = Gridap.Visualization.VisualizationGrid(n_trian_tag, n_ref_grids0)
    n_Γ_data = Gridap.Visualization._prepare_pdata(n_trian_tag, n_Γ_dict, n_visgrid0.cell_to_refpoints)
    airfoil_normals = n_Γ_data["n_Γ"]
    

    airfoil_normals_top = airfoil_normals[IDX_TOP]
    airfoil_normals_bottom = airfoil_normals[IDX_BOTTOM]
    return airfoil_normals_top,airfoil_normals_bottom
end



function ParametersAdj.AirfoilModel(model, mcase::Airfoil;  tag="airfoil")
    @sunpack AoA = mcase
    nodesu,nodesl,params = get_nodes_idx(model, AoA,tag)

    ap = AirfoilPoints(getindex.(nodesu,1),getindex.(nodesl,1),getindex.(nodesu,2),getindex.(nodesl,2))
    nu,nl = get_normals(model, params, tag)
    an = AirfoilNormals(map(n-> [n...], nu),map(n-> [n...], nl))

    return AirfoilModel(ap,model,an,params)
end


function get_aerodynamic_features(am::AirfoilModel, uh,ph;tag="airfoil")

    @unpack IDX_TOP_UNIQUE,IDX_BOTTOM_UNIQUE=am.params
    u_in = 1.0
    q = 0.5 .* u_in^2

    f = (reffe) -> Gridap.Geometry.UnstructuredGrid(reffe)

    Γ = BoundaryTriangulation(am.model; tags=tag)
    cf = Dict("uh"=>uh,"ph"=>ph)
    trian_tag = Γ.trian
    
    ref_grids0 = map(f, Gridap.Geometry.get_reffes(trian_tag))
    visgrid0 = Gridap.Visualization.VisualizationGrid(trian_tag, ref_grids0)
    

    
    pdata = Gridap.Visualization._prepare_pdata(Γ,cf,visgrid0.cell_to_refpoints)

    
    pressure_top = pdata["ph"][IDX_TOP_UNIQUE]
    pressure_bottom  = pdata["ph"][IDX_BOTTOM_UNIQUE]
    cp_top = pressure_top ./ q
    cp_bottom = pressure_bottom./ q

     return AirfoilScalar(cp_top,cp_bottom)
end

