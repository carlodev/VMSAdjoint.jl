#############################################################################
# Mesh deformation (remeshing-free shape perturbation)
#############################################################################
#
# The shape-sensitivity loop needs, for every design parameter, the mesh
# corresponding to a perturbed geometry. Regenerating it with Gmsh costs a full
# mesh generation per design parameter per optimizer iteration, is serial, and
# — for unstructured meshes — does not preserve the mesh topology, which is why
# CST designs were previously restricted to structured meshes.
#
# Deforming the baseline mesh instead removes all three problems:
#
#   1. no Gmsh call in the inner loop;
#   2. identical topology by construction, so the DOF-wise mesh velocity
#      (m1 - m0)/δ is exact rather than an artifact of two independent meshings;
#   3. the operator of the mesh-motion problem does not depend on the design
#      perturbation — only its Dirichlet data does — so it is assembled and
#      factorized ONCE and reused for all `Ndes` perturbations.
#
# The extension is a stiffened harmonic (Laplace) problem
#
#     ∇·(μ ∇d) = 0    in Ω,     d = d_surface on Γ_airfoil,   d = 0 on the farfield
#
# with cell-wise stiffening μ = (V_max/V_cell)^χ, which makes the small
# boundary-layer cells nearly rigid so the near-wall mesh translates with the
# surface instead of absorbing all the strain. This is the standard
# Jacobian-stiffened pseudo-elastic mesh motion, in its cheapest (scalar
# Laplacian) form.
#############################################################################

"""
    chordwise_interpolant(x, y)

Linear interpolant of `y(x)` over airfoil surface stations, tolerant of how the
station arrays actually arrive.

`AirfoilPoints` are not guaranteed to be sorted or unique: coordinate files are
commonly ordered trailing-edge → leading-edge → trailing-edge (so one side runs
backwards), and the leading/trailing edge point is frequently repeated. The
stations are therefore sorted and de-duplicated before the interpolant is built.
Dropping a repeated station is safe here: both copies sit at the same chordwise
location, where the surface displacement of either parametrization is pinned.

Extrapolation is flat, so query points marginally outside the sampled chord
(mesh nodes at a rounded leading edge) pick up the nearest station's value
rather than diverging.
"""
function chordwise_interpolant(x::AbstractVector, y::AbstractVector)
    length(x) == length(y) || throw(ArgumentError(
        "chordwise stations and values have different lengths ($(length(x)) vs $(length(y)))"))

    perm = sortperm(x)
    xs, ys = collect(float.(x[perm])), collect(float.(y[perm]))

    keep = [1]
    for i in 2:length(xs)
        xs[i] > xs[keep[end]] + 1e-12 && push!(keep, i)
    end
    length(keep) >= 2 || throw(ArgumentError(
        "airfoil surface has fewer than 2 distinct chordwise stations"))

    return linear_interpolation(xs[keep], ys[keep], extrapolation_bc=Flat())
end

"""
    airfoil_surface_displacement(ap0::AirfoilPoints, ap1::AirfoilPoints, AoA)

Displacement field on the airfoil surface taking geometry `ap0` to `ap1`, as a
callable `x -> VectorValue{2}` in the *mesh* frame.

Both supported parametrizations perturb a design parameter by a purely *vertical*
displacement of the surface in the airfoil-fixed design frame: CST changes a
weight, which moves `y(x)` at fixed `x`; RBF morphs with an explicitly vertical
`[0, δ]` displacement. So the surface displacement is fully described by `Δy(x)`
per side, and is rotated into the mesh frame by the angle of attack.

The two geometries need not share a chordwise distribution — `Δy` is formed by
resampling `ap1` onto the baseline stations of `ap0` — which matters for CST,
whose `perturb_DesignParameter` rebuilds the point set from the weights rather
than reusing the baseline `x` array.

The upper/lower side of a query point is decided by which baseline surface it is
closer to, not by the sign of `y`: on a cambered airfoil the lower surface can
sit above `y = 0`.
"""
function airfoil_surface_displacement(ap0::AirfoilPoints, ap1::AirfoilPoints, AoA::Real)
    itp = chordwise_interpolant

    yu1, yl1 = itp(ap1.xu, ap1.yu), itp(ap1.xl, ap1.yl)
    Δyu = itp(ap0.xu, yu1.(ap0.xu) .- ap0.yu)
    Δyl = itp(ap0.xl, yl1.(ap0.xl) .- ap0.yl)
    yu0, yl0 = itp(ap0.xu, ap0.yu), itp(ap0.xl, ap0.yl)

    AoAf = Float64(AoA)

    return function (x)
        xd, yd = rotate_points([x[1], x[2]], -AoAf)          # mesh frame -> design frame
        Δy = abs(yd - yu0(xd)) <= abs(yd - yl0(xd)) ? Δyu(xd) : Δyl(xd)
        dx, dy = rotate_points([0.0, Float64(Δy)], AoAf)     # design frame -> mesh frame
        return VectorValue(dx, dy)
    end
end

#############################################################################
# Mesh-motion operator (assembled and factorized once per baseline mesh)
#############################################################################

"""
    MeshMotionCache

Factorized stiffened-Laplacian mesh-motion operator for one baseline mesh, plus
the spaces needed to apply it. Built by [`mesh_motion_cache`](@ref) and reused
for every design perturbation of that mesh — the matrix is independent of the
perturbation, only the Dirichlet data changes.

- `Vd`       — test space, Dirichlet on airfoil + farfield
- `VV1`      — unconstrained order-1 vector space; its free DOFs are the mesh
               nodes, in node order (see [`node_ordered_dofs_check`](@ref))
- `a`        — the stiffened bilinear form
- `ns`       — numerical setup (factorization) of the free-free block
- `n_nodes`  — number of mesh nodes
"""
struct MeshMotionCache{TV,TW,TA,TN}
    Vd::TV
    VV1::TW
    a::TA
    ns::TN
    n_nodes::Int
end

"""
    node_ordered_dofs_check(model, VV1, D)

Assert that the free DOFs of the unconstrained order-1 vector Lagrangian space
`VV1` are laid out as `dof = (node-1)*D + component`, i.e. that
`reshape(get_free_dof_values(...), D, :)` is exactly the mesh node-coordinate
array.

This layout is what lets the mesh-motion solution be written straight onto the
mesh nodes with no permutation. It holds for Gridap's Lagrangian spaces, but it
is an internal convention, so it is verified rather than assumed — the check is
one pass over the nodes and runs once per baseline mesh.
"""
function node_ordered_dofs_check(model, VV1, D::Int)
    node_coords = Gridap.Geometry.get_node_coordinates(Gridap.Geometry.get_grid(model))
    m0 = get_free_dof_values(interpolate_everywhere(x -> x, VV1))
    length(m0) == D * length(node_coords) || error(
        "Mesh-motion space has $(length(m0)) DOFs for $(length(node_coords)) nodes in $D dimensions; " *
        "expected a node-major order-1 Lagrangian layout.")
    M = reshape(m0, D, :)
    for (n, xc) in enumerate(node_coords)
        for c in 1:D
            isapprox(M[c, n], xc[c]; atol=1e-10) || error(
                "Gridap DOF ordering is not node-major for this mesh (node $n, component $c): " *
                "the mesh-motion node update would scatter displacements to the wrong nodes.")
        end
    end
    return length(node_coords)
end

"""
    stiffening_field(Ω, χ)

Cell-wise mesh stiffness `μ = (V_max / V_cell)^χ`: small cells (the boundary
layer) are stiff and move rigidly with the surface, large far-field cells absorb
the strain. `χ = 0` gives an unstiffened harmonic extension.
"""
function stiffening_field(Ω, χ::Float64)
    vols = get_cell_measure(Ω)
    vmax = maximum(vols)
    return CellField([(vmax / v)^χ for v in vols], Ω)
end

"""
    mesh_motion_cache(model, D, χ=1.0; solver=LUSolver())

Assemble and factorize the mesh-motion operator of `model`. The result is
independent of any particular design perturbation and is meant to be built once
per baseline mesh and cached in `am.params[:mesh_motion]` (see
[`mesh_motion_cache!`](@ref)).
"""
function mesh_motion_cache(model, D::Int, χ::Float64=1.0; solver=LUSolver())
    reffe1 = ReferenceFE(lagrangian, VectorValue{D,Float64}, 1)

    Vd = TestFESpace(model, reffe1, conformity=:H1,
        dirichlet_tags=["airfoil", "inlet", "outlet", "limits"])
    VV1 = FESpace(model, reffe1; conformity=:H1)
    n_nodes = node_ordered_dofs_check(model, VV1, D)

    Ω = Triangulation(model)
    dΩ = Measure(Ω, 2)
    μ = stiffening_field(Ω, χ)
    a(d, w) = ∫(μ * (∇(d) ⊙ ∇(w)))dΩ

    # homogeneous trial space: the Dirichlet lifting is applied per perturbation
    # in `solve_mesh_displacement`, so this matrix (and its factorization) is
    # shared by every design parameter
    zero_disp = VectorValue(zeros(D)...)
    U0 = TrialFESpace(Vd, [zero_disp, zero_disp, zero_disp, zero_disp])
    A = assemble_matrix(a, U0, Vd)
    ns = numerical_setup(symbolic_setup(solver, A), A)

    return MeshMotionCache(Vd, VV1, a, ns, n_nodes)
end

"""
    mesh_motion_cache!(am::AirfoilModel, D, χ=1.0)

Return the [`MeshMotionCache`](@ref) of `am`, building and registering it in
`am.params[:mesh_motion]` on first use.
"""
function mesh_motion_cache!(am::AirfoilModel, D::Int, χ::Float64=1.0)
    if !haskey(am.params, :mesh_motion)
        updatekey(am.params, :mesh_motion, mesh_motion_cache(am.model, D, χ))
    end
    return am.params[:mesh_motion]
end

"""
    solve_mesh_displacement(cache::MeshMotionCache, dispfun, D)

Extend the airfoil-surface displacement `dispfun` (a callable `x -> VectorValue`)
into the volume by solving the cached stiffened-Laplacian problem. Returns the
displacement as an `FEFunction` on the unconstrained space `cache.VV1`, so it
carries the Dirichlet values too and can be evaluated anywhere in Ω (including
on the airfoil boundary, where the sensitivity kernel needs it).

Only the right-hand side is assembled here; the factorization is reused.
"""
function solve_mesh_displacement(cache::MeshMotionCache, dispfun, D::Int)
    @unpack Vd, VV1, a, ns = cache

    zero_disp = VectorValue(zeros(D)...)
    Ug = TrialFESpace(Vd, [dispfun, zero_disp, zero_disp, zero_disp])

    # dg carries the Dirichlet data with zero free values: A·d_free = -a(dg, ·).
    # Gridap has no unary minus on a DomainContribution, so the lifting is
    # assembled first and negated on the vector.
    dg = FEFunction(Ug, zeros(num_free_dofs(Ug)))
    b = assemble_vector(w -> a(dg, w), Vd)
    rmul!(b, -1)

    d_free = similar(b)
    fill!(d_free, zero(eltype(d_free)))
    solve!(d_free, ns, b)

    dh = FEFunction(Ug, d_free)
    return interpolate_everywhere(dh, VV1)
end

#############################################################################
# Applying the displacement to the mesh
#############################################################################

"""
    cell_signed_areas(model)

Signed area of every cell of a 2D first-order mesh, from the shoelace formula on
its vertices. `nothing` for meshes this cannot handle (non-2D, or cell types
other than TRI/QUAD), in which case [`deform_model`](@ref) falls back to the
volume-ratio heuristic alone.

Used to detect mesh tangling exactly: a valid deformed mesh has all signed areas
of the same sign as the baseline.
"""
function cell_signed_areas(model)
    grid = Gridap.Geometry.get_grid(model)
    num_cell_dims(grid) == 2 && num_point_dims(grid) == 2 || return nothing

    coords = Gridap.Geometry.get_node_coordinates(grid)
    cell_nodes = Gridap.Geometry.get_cell_node_ids(grid)

    areas = Vector{Float64}(undef, length(cell_nodes))
    for (c, nodes) in enumerate(cell_nodes)
        # Gridap orders QUAD nodes tensor-product-wise (00,10,01,11); the
        # shoelace formula needs them cyclic
        ring = length(nodes) == 3 ? nodes :
               length(nodes) == 4 ? nodes[[1, 2, 4, 3]] : return nothing
        s = 0.0
        for k in eachindex(ring)
            p = coords[ring[k]]
            q = coords[ring[mod1(k + 1, length(ring))]]
            s += p[1] * q[2] - q[1] * p[2]
        end
        areas[c] = 0.5 * s
    end
    return areas
end

"""
    deform_model(model, dh, D; min_volume_ratio=0.1, check_tangling=true)

Displace the nodes of `model` by the order-1 displacement field `dh` (as
returned by [`solve_mesh_displacement`](@ref)) and return the deformed
`DiscreteModel`. Topology, cell types and the face labeling are shared with the
baseline, so the deformed model has exactly the same DOF layout — which is what
makes the DOF-wise mesh velocity `(m1 - m0)/δ` meaningful.

Mesh validity is checked before the deformed model is returned:

- 2D TRI/QUAD meshes: every cell must keep the sign of its signed area
  (exact tangling detection);
- always: no cell may shrink below `min_volume_ratio` of its baseline volume,
  nor grow beyond `1/min_volume_ratio`.

An invalid deformation raises rather than silently producing a garbage
sensitivity. If a legitimate design step trips it, the displacement is too large
for a single extension — reduce the step or fall back to remeshing.
"""
function deform_model(model, dh, D::Int; min_volume_ratio::Float64=0.1, check_tangling::Bool=true)
    grid = Gridap.Geometry.get_grid(model)
    node_coords = Gridap.Geometry.get_node_coordinates(grid)

    disp = reshape(get_free_dof_values(dh), D, :)
    size(disp, 2) == length(node_coords) || error(
        "Displacement field has $(size(disp,2)) nodes, mesh has $(length(node_coords))")

    new_coords = [node_coords[n] + VectorValue(ntuple(c -> disp[c, n], D))
                  for n in eachindex(node_coords)]

    new_grid = Gridap.Geometry.UnstructuredGrid(
        new_coords,
        Gridap.Geometry.get_cell_node_ids(grid),
        Gridap.Geometry.get_reffes(grid),
        Gridap.Geometry.get_cell_type(grid),
        Gridap.Geometry.OrientationStyle(grid))

    topo = Gridap.Geometry.get_grid_topology(model)
    new_topo = Gridap.Geometry.UnstructuredGridTopology(
        new_grid,
        Gridap.Geometry.get_faces(topo, num_cell_dims(topo), 0),
        Gridap.Geometry.get_vertex_node(model))

    new_model = Gridap.Geometry.UnstructuredDiscreteModel(
        new_grid, new_topo, get_face_labeling(model))

    validate_deformed_model(model, new_model, min_volume_ratio, check_tangling)

    return new_model
end

"Raise if the deformed mesh is tangled or severely distorted relative to the baseline."
function validate_deformed_model(model0, model1, min_volume_ratio::Float64, check_tangling::Bool)
    if check_tangling
        a0 = cell_signed_areas(model0)
        a1 = cell_signed_areas(model1)
        if !isnothing(a0) && !isnothing(a1)
            tangled = findall(i -> sign(a1[i]) != sign(a0[i]), eachindex(a0))
            isempty(tangled) || error(
                "Mesh deformation tangled $(length(tangled)) cell(s) (first: cell $(first(tangled))): " *
                "the surface displacement is too large for a single harmonic extension.")
        end
    end

    v0 = get_cell_measure(Triangulation(model0))
    v1 = get_cell_measure(Triangulation(model1))
    ratios = v1 ./ v0
    lo, hi = extrema(ratios)
    (lo >= min_volume_ratio && hi <= 1 / min_volume_ratio) || error(
        "Mesh deformation distorted cell volumes by a factor in [$(round(lo, digits=4)), " *
        "$(round(hi, digits=4))], outside the accepted " *
        "[$min_volume_ratio, $(1/min_volume_ratio)]: reduce the design step or remesh.")

    return nothing
end

#############################################################################
# Entry point
#############################################################################

"""
    deform_airfoil_mesh(am0, ap0, ap1, simcase; χ=1.0, min_volume_ratio=0.1)

Deform the mesh of `am0` so its airfoil boundary follows the geometry change
from `ap0` to `ap1`. Returns `(model1, dh)`:

- `model1` — the deformed `DiscreteModel`, topologically identical to `am0.model`
- `dh`     — the volume displacement field, an `FEFunction` valid over all of Ω
             (including the airfoil boundary), which the sensitivity kernel
             consumes directly as the mesh-velocity field.

`ap0`/`ap1` must be the *design* point sets (`adesign.ap`), expressed in the
airfoil-fixed design frame on a common chordwise distribution — NOT
`am0.ap`, which holds the mesh surface nodes already rotated into the mesh
frame by the angle of attack.

The mesh-motion operator is cached on `am0` after the first call, so a full
`Ndes`-parameter sweep pays one assembly and one factorization in total.
"""
function deform_airfoil_mesh(am0::AirfoilModel, ap0::AirfoilPoints, ap1::AirfoilPoints,
    simcase::Airfoil; χ::Float64=1.0, min_volume_ratio::Float64=0.1)
    @sunpack D, AoA = simcase
    D == 2 || throw(ArgumentError(
        "Mesh deformation currently builds the surface displacement from 2D airfoil " *
        "sections; 3D needs a spanwise-aware surface parametrization."))

    cache = mesh_motion_cache!(am0, D, χ)
    dispfun = airfoil_surface_displacement(ap0, ap1, AoA)
    dh = solve_mesh_displacement(cache, dispfun, D)
    model1 = deform_model(am0.model, dh, D; min_volume_ratio=min_volume_ratio)

    return model1, dh
end
