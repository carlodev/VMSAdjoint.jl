# Unit tests for the remeshing-free shape-perturbation path.
#
# These run on a synthetic square domain whose four boundary edges are relabelled
# with the airfoil-case tag names, so the mesh-motion operator can be exercised
# without generating an airfoil mesh.

using Gridap
using Gridap.Geometry
using Gridap.FESpaces
using AirfoilTools
using Gridap.Geometry: get_node_coordinates, num_nodes, get_cell_node_ids, get_grid
using VMSAdjoint.Interfaces: airfoil_surface_displacement, mesh_motion_cache,
    solve_mesh_displacement, deform_model, cell_signed_areas,
    node_ordered_dofs_check, stiffening_field, chordwise_interpolant

@testset "chordwise_interpolant" begin
    x = collect(LinRange(0.0, 1.0, 5))
    y = 2 .* x

    # airfoil coordinate files commonly run trailing-edge -> leading-edge, i.e.
    # one surface arrives with decreasing x
    f_rev = chordwise_interpolant(reverse(x), reverse(y))
    @test f_rev(0.25) ≈ 0.5
    @test f_rev(0.6) ≈ 1.2

    # a repeated station (blunt/shared leading or trailing edge) must not throw
    f_dup = chordwise_interpolant(vcat(x, 1.0), vcat(y, 2.0))
    @test f_dup(0.5) ≈ 1.0

    # flat extrapolation outside the sampled chord
    @test f_rev(-0.3) ≈ 0.0
    @test f_rev(1.7) ≈ 2.0

    @test_throws ArgumentError chordwise_interpolant([0.0, 1.0], [0.0])
    @test_throws ArgumentError chordwise_interpolant([0.5, 0.5], [1.0, 2.0])
end

"Square model with the airfoil-case boundary tags mapped onto its four edges."
function tagged_square(n::Int=6)
    model = simplexify(CartesianDiscreteModel((0, 1, 0, 1), (n, n)))
    labels = get_face_labeling(model)
    # 2D Cartesian boundary entities: tags 1-4 are corners, 5-8 the edges
    # (x-min, x-max, y-min, y-max)
    add_tag_from_tags!(labels, "airfoil", ["tag_6"])
    add_tag_from_tags!(labels, "inlet", ["tag_5"])
    add_tag_from_tags!(labels, "outlet", ["tag_7"])
    add_tag_from_tags!(labels, "limits", ["tag_8"])
    return model
end

@testset "airfoil_surface_displacement" begin
    x = collect(LinRange(0.0, 1.0, 11))
    ap0 = AirfoilPoints(x, x, 0.1 .* x, -0.1 .* x)
    # upper surface raised by 0.02, lower surface untouched
    ap1 = AirfoilPoints(x, x, 0.1 .* x .+ 0.02, -0.1 .* x)

    f = airfoil_surface_displacement(ap0, ap1, 0.0)
    # a point sitting on the upper surface picks up the upper-side Δy
    @test f(VectorValue(0.5, 0.05))[2] ≈ 0.02
    @test f(VectorValue(0.5, 0.05))[1] ≈ 0.0
    # a point on the lower surface is unmoved
    @test f(VectorValue(0.5, -0.05))[2] ≈ 0.0

    # side selection is by proximity to the baseline surfaces, not by sign(y):
    # a cambered geometry whose lower surface sits entirely above y = 0
    ap0c = AirfoilPoints(x, x, 0.3 .+ 0.1 .* x, 0.1 .+ 0.1 .* x)
    ap1c = AirfoilPoints(x, x, 0.3 .+ 0.1 .* x, 0.1 .+ 0.1 .* x .- 0.05)
    fc = airfoil_surface_displacement(ap0c, ap1c, 0.0)
    @test fc(VectorValue(0.5, 0.15))[2] ≈ -0.05   # on the lower surface, y > 0
    @test fc(VectorValue(0.5, 0.35))[2] ≈ 0.0     # on the upper surface

    # the pure-y design perturbation is rotated into the mesh frame
    frot = airfoil_surface_displacement(ap0, ap1, 90.0)
    d = frot(rotate_points([0.5, 0.05], 90.0) |> v -> VectorValue(v[1], v[2]))
    @test d[1] ≈ 0.02 atol = 1e-12
    @test d[2] ≈ 0.0 atol = 1e-12

    # the perturbed geometry may carry its own chordwise distribution (CST
    # regenerates one from the weights): Δy is resampled onto the baseline
    # stations, so a denser point set describing the SAME shape gives no motion
    xd = collect(LinRange(0.0, 1.0, 37))
    ap1_dense = AirfoilPoints(xd, xd, 0.1 .* xd .+ 0.02, -0.1 .* xd)
    fd = airfoil_surface_displacement(ap0, ap1_dense, 0.0)
    @test fd(VectorValue(0.5, 0.05))[2] ≈ 0.02
    ap_same_dense = AirfoilPoints(xd, xd, 0.1 .* xd, -0.1 .* xd)
    fsame = airfoil_surface_displacement(ap0, ap_same_dense, 0.0)
    @test abs(fsame(VectorValue(0.5, 0.05))[2]) < 1e-12
    @test abs(fsame(VectorValue(0.3, -0.03))[2]) < 1e-12
end

@testset "node/DOF layout" begin
    model = tagged_square()
    VV1 = FESpace(model, ReferenceFE(lagrangian, VectorValue{2,Float64}, 1); conformity=:H1)
    # the check returns the node count when the layout is node-major
    @test node_ordered_dofs_check(model, VV1, 2) == num_nodes(model)
end

@testset "stiffening_field" begin
    model = tagged_square()
    Ω = Triangulation(model)
    dΩ = Measure(Ω, 2)
    # uniform mesh -> uniform stiffness = 1 everywhere
    μ = stiffening_field(Ω, 1.0)
    @test sum(∫(μ)dΩ) ≈ 1.0 rtol = 1e-10
    # χ = 0 disables stiffening
    @test sum(∫(stiffening_field(Ω, 0.0))dΩ) ≈ 1.0 rtol = 1e-10
end

@testset "mesh motion solve" begin
    model = tagged_square(8)
    cache = mesh_motion_cache(model, 2)

    disp(x) = VectorValue(0.03 * sin(pi * x[2]), 0.0)
    dh = solve_mesh_displacement(cache, disp, 2)

    Γa = BoundaryTriangulation(model; tags="airfoil")
    Γi = BoundaryTriangulation(model; tags="inlet")
    dΓa, dΓi = Measure(Γa, 4), Measure(Γi, 4)

    # Dirichlet data is reproduced on the airfoil and homogeneous on the farfield
    err = dh - CellField(disp, Γa)
    @test sqrt(abs(sum(∫(err ⋅ err)dΓa))) < 1e-10
    @test sqrt(abs(sum(∫(dh ⋅ dh)dΓi))) < 1e-10

    # a harmonic extension cannot exceed its boundary data (maximum principle)
    vals = reshape(get_free_dof_values(dh), 2, :)
    @test maximum(abs, vals) <= 0.03 + 1e-10
end

@testset "deform_model" begin
    model = tagged_square(6)
    grid = get_grid(model)
    coords0 = copy(get_node_coordinates(grid))

    VV1 = FESpace(model, ReferenceFE(lagrangian, VectorValue{2,Float64}, 1); conformity=:H1)

    # a prescribed nodal field must land on exactly the nodes it names: this is
    # the end-to-end check of the node-major DOF assumption
    shift = [VectorValue(0.001 * n, -0.002 * n) for n in 1:length(coords0)]
    free = zeros(2 * length(coords0))
    for n in eachindex(shift)
        free[2*(n-1)+1] = shift[n][1]
        free[2*(n-1)+2] = shift[n][2]
    end
    dh = FEFunction(VV1, free)

    model1 = deform_model(model, dh, 2)
    coords1 = get_node_coordinates(get_grid(model1))
    @test all(coords1[n] ≈ coords0[n] + shift[n] for n in eachindex(coords0))

    # topology, cell types and labels are shared with the baseline
    @test get_cell_node_ids(get_grid(model1)) == get_cell_node_ids(grid)
    @test num_cells(model1) == num_cells(model)
    @test num_nodes(model1) == num_nodes(model)

    # a rigid translation leaves every cell area unchanged
    dh_rigid = FEFunction(VV1, repeat([0.01, 0.02], length(coords0)))
    a0 = cell_signed_areas(model)
    a1 = cell_signed_areas(deform_model(model, dh_rigid, 2))
    @test a0 ≈ a1
end

@testset "deform_model rejects invalid meshes" begin
    model = tagged_square(4)
    VV1 = FESpace(model, ReferenceFE(lagrangian, VectorValue{2,Float64}, 1); conformity=:H1)
    coords0 = get_node_coordinates(get_grid(model))

    # collapse everything onto a point: cells are destroyed, must not pass
    free = zeros(2 * length(coords0))
    for n in eachindex(coords0)
        free[2*(n-1)+1] = -coords0[n][1]
        free[2*(n-1)+2] = -coords0[n][2]
    end
    @test_throws ErrorException deform_model(model, FEFunction(VV1, free), 2)

    # mirror the mesh through x -> -x: every cell inverts, tangling detection
    # must catch it even though the areas keep their magnitude
    free_flip = zeros(2 * length(coords0))
    for n in eachindex(coords0)
        free_flip[2*(n-1)+1] = -2 * coords0[n][1]
    end
    @test_throws ErrorException deform_model(model, FEFunction(VV1, free_flip), 2)
end
