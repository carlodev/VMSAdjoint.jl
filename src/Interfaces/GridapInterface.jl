using AirfoilTools
using Gridap
using Parameters
using SegregatedVMSSolver

#############################################################################
# Boundary point extraction helpers
#############################################################################

"Visualization grid of a (boundary) triangulation, used to sample fields at its nodes."
function _visualization_grid(trian)
    ref_grids = map(reffe -> Gridap.Geometry.UnstructuredGrid(reffe), Gridap.Geometry.get_reffes(trian))
    return Gridap.Visualization.VisualizationGrid(trian, ref_grids)
end

"Sample `fields` (a `Dict{String,<:Any}`) at the visualization nodes of `trian`."
function _point_data(trian, fields::Dict, visgrid=_visualization_grid(trian))
    return Gridap.Visualization._prepare_pdata(trian, fields, visgrid.cell_to_refpoints)
end

"""
    boundary_outward_normals(model, tag)

Outward-pointing normals of the boundary `tag`, sampled at its visualization nodes.
"""
function boundary_outward_normals(model, tag::String)
    Γ = BoundaryTriangulation(model; tags=tag)
    n_Γ = -1 .* get_normal_vector(Γ) # sign flip: Gridap normals point into the body
    return _point_data(n_Γ.trian, Dict("n_Γ" => n_Γ))["n_Γ"]
end

"""
    airfoil_boundary!(am::AirfoilModel, order::Int)

Return `(dΓ, nΓ)` of the airfoil boundary, with `nΓ` pointing outward from the
body, registering `:Γ/:dΓ/:nΓ` in `am.params` on first use. Single shared
construction for the force computation, the sensitivity kernel and the solvers.
"""
function airfoil_boundary!(am::AirfoilModel, order::Int)
    params = am.params
    if !haskey(params, :dΓ)
        Γ = BoundaryTriangulation(am.model; tags="airfoil")
        updatekey(params, :Γ, Γ)
        updatekey(params, :dΓ, Measure(Γ, 2 * order))
        updatekey(params, :nΓ, -get_normal_vector(Γ)) # sign flip: Gridap normals point into the body
    end
    return params[:dΓ], params[:nΓ]
end

#############################################################################
# Airfoil point/normal identification
#############################################################################

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

"""
    remove_outliers_from_curve(points, tol)

Drop interior points deviating more than `tol` from the midpoint of their
neighbours; endpoints are always kept. Returns `(filtered_points, kept_indices)`.
"""
function remove_outliers_from_curve(points::Vector, tol::Float64)
    filtered_pts = [points[1]] # keep first point
    idxs = [1]

    for i in 2:length(points)-1
        expected = 0.5 .* (points[i-1] + points[i+1])
        if norm(points[i] - expected) < tol
            push!(filtered_pts, points[i])
            push!(idxs, i)
        end
    end

    push!(filtered_pts, points[end]) # keep last point
    push!(idxs, length(points))

    return filtered_pts, idxs
end

"""
    get_nodes_idx(model, AoA, tag)

Identify and sort the airfoil surface nodes of `model`, splitting them into
top/bottom sides by the sign of the outward normal. Returns the top and bottom
point sets plus the index bookkeeping needed to sample fields on them.
"""
function get_nodes_idx(model, AoA::Float64, tag::String)
    Γ = BoundaryTriangulation(model; tags=tag)
    visgrid = _visualization_grid(Γ.trian)
    airfoil_points0 = visgrid.sub_grid.node_coordinates
    idx_uniques = uniqueidx(airfoil_points0)

    normals = boundary_outward_normals(model, tag)

    airfoil_points = airfoil_points0[idx_uniques]
    airfoil_normals = normals[idx_uniques]

    idx_top = findall(n -> n[2] > 0.1, airfoil_normals)
    idx_bottom = findall(n -> n[2] < -0.1, airfoil_normals)

    IDX_TOP = idx_top[sortperm(getindex.(airfoil_points[idx_top], 1))]
    IDX_BOTTOM = idx_bottom[sortperm(getindex.(airfoil_points[idx_bottom], 1))]
    IDX_TOP_UNIQUE = idx_uniques[IDX_TOP]
    IDX_BOTTOM_UNIQUE = idx_uniques[IDX_BOTTOM]

    airfoil_points_top, idx_top = remove_outliers_from_curve(airfoil_points[IDX_TOP], 0.01)
    airfoil_points_bottom, idx_bottom = remove_outliers_from_curve(airfoil_points[IDX_BOTTOM], 0.01)

    params = Dict(
        :IDX_TOP_UNIQUE => IDX_TOP_UNIQUE[idx_top],
        :IDX_BOTTOM_UNIQUE => IDX_BOTTOM_UNIQUE[idx_bottom],
        :IDX_TOP => IDX_TOP, :IDX_BOTTOM => IDX_BOTTOM)

    return airfoil_points_top, airfoil_points_bottom, params
end

"Airfoil outward normals split into top/bottom, in the ordering of `get_nodes_idx`."
function get_normals(model, params, tag)
    @unpack IDX_TOP, IDX_BOTTOM = params
    normals = boundary_outward_normals(model, tag)
    return normals[IDX_TOP], normals[IDX_BOTTOM]
end

function ParametersAdj.AirfoilModel(model, mcase::Airfoil; tag="airfoil")
    @sunpack AoA, D = mcase
    D == 2 || throw(ArgumentError(
        "3D meshes are not yet supported: mesh generation and airfoil surface " *
        "extraction are 2D-only (the solver and sensitivity kernels are dimension-generic)."))
    nodesu, nodesl, params = get_nodes_idx(model, AoA, tag)

    ap = AirfoilPoints(getindex.(nodesu, 1), getindex.(nodesl, 1), getindex.(nodesu, 2), getindex.(nodesl, 2))
    nu, nl = get_normals(model, params, tag)
    an = AirfoilNormals(map(n -> [n...], nu), map(n -> [n...], nl))

    return AirfoilModel(ap, model, an, params)
end

#############################################################################
# Aerodynamic post-processing
#############################################################################

"""
    get_aerodynamic_features(am::AirfoilModel, physicalp::PhysicalParameters, uh, ph; tag="airfoil")

Pressure coefficient distribution on the top/bottom airfoil surfaces,
normalized by the dynamic-pressure reference of `physicalp` (consistent with
`compute_airfoil_coefficients`).
"""
function get_aerodynamic_features(am::AirfoilModel, physicalp::PhysicalParameters, uh, ph; tag="airfoil")
    @unpack IDX_TOP_UNIQUE, IDX_BOTTOM_UNIQUE = am.params
    @unpack c, u_in_mag = physicalp
    q = 0.5 * c * u_in_mag^2

    Γ = BoundaryTriangulation(am.model; tags=tag)
    visgrid = _visualization_grid(Γ.trian)
    pdata = _point_data(Γ, Dict("uh" => uh, "ph" => ph), visgrid)

    cp_top = pdata["ph"][IDX_TOP_UNIQUE] ./ q
    cp_bottom = pdata["ph"][IDX_BOTTOM_UNIQUE] ./ q

    return AirfoilScalar(cp_top, cp_bottom)
end

"""
    get_aerodynamic_features(am::AirfoilModel, physicalp::PhysicalParameters, uh; tag="airfoil")

Skin-friction coefficient distribution on the top/bottom airfoil surfaces
(the pressure is not needed, hence the shorter argument list: `uh, ph` gives Cp,
`uh` alone gives Cf).

The wall shear is the tangential part of the viscous traction,

    τ_w = (I - n ⊗ n) ⋅ [ν (∇u + ∇uᵀ) ⋅ n],

i.e. the viscous term of [`compute_airfoil_forces`](@ref) with its normal
component removed, normalized by the same dynamic-pressure reference as the Cp
method. The sign comes from the projection of `τ_w` on the freestream
direction, so `cf < 0` marks reversed (separated) flow. Like
[`tangential_normal_derivative`](@ref), the projector form carries over to 3D,
where `|τ_w|` also picks up the spanwise shear.
"""
function get_aerodynamic_features(am::AirfoilModel, physicalp::PhysicalParameters, uh; tag="airfoil")
    @unpack IDX_TOP_UNIQUE, IDX_BOTTOM_UNIQUE = am.params
    @unpack c, u_in, u_in_mag, ν = physicalp
    q = 0.5 * c * u_in_mag^2

    Γ = BoundaryTriangulation(am.model; tags=tag)
    nΓ = -get_normal_vector(Γ) # sign flip: Gridap normals point into the body
    tν = ν * (∇(uh) + transpose(∇(uh))) ⋅ nΓ # viscous traction, as in compute_airfoil_forces
    τw = tν - (tν ⋅ nΓ) * nΓ

    visgrid = _visualization_grid(Γ.trian)
    pdata = _point_data(Γ, Dict("τw" => τw), visgrid)

    D = num_cell_dims(am.model)
    ê = VectorValue(ntuple(i -> i <= length(u_in) ? u_in[i] : 0.0, D)) / u_in_mag
    cf(τ) = sign(τ ⋅ ê) * norm(τ) / q

    cf_top = cf.(pdata["τw"][IDX_TOP_UNIQUE])
    cf_bottom = cf.(pdata["τw"][IDX_BOTTOM_UNIQUE])

    return AirfoilScalar(cf_top, cf_bottom)
end
