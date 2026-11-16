


"""
    set_curvature_sizing!(θ_max_deg, h_min)

Limit every mesh edge on a curved boundary to a turning angle of at most
`θ_max_deg` degrees, with `h_min` as the absolute floor. `θ_max_deg = 0`
disables it and reproduces uniform sizing exactly.

Gmsh expresses this as "elements per 2π of curvature", hence `360/θ_max_deg`.

Why this exists
---------------
The airfoil *spline* is well resolved — `n0012.csv` is cosine-clustered near the
nose (x = 0, 5.8e-4, 2.3e-3, 5.3e-3, …). The *mesh* was not. The 129 CSV points
are `addSpline` control points, so they bound no curve of their own and Gmsh
ignores the per-point size given to `addPoint`; the whole airfoil is one NURBS
whose sizing came only from its two endpoint vertices, both at `rpoint`. The
surface therefore came out uniform at 0.005c: on a NACA0012 (r_LE = 0.0159c)
only ~7 boundary nodes landed inside one leading-edge radius, each nose facet
turning through ~18°, against <0.3° at mid-chord.

That matters because the boundary is a polygon of straight-sided triangles, so
`nΓ` is piecewise constant, and ~18° per facet puts an O(10°) error on the
normal exactly at the nose. `nΓ` enters the shape-sensitivity kernel three
times, and once first-order and non-cancelling, as `δβ = (dh⋅nΓ)/δ`
([`deformation_normal_field`](@ref)). The higher-derivative-free kernel of
[`compute_gradient`](@ref) is exact in the continuum — its pressure and
curvature terms cancel analytically, using `∇u = (∂ₙu)⊗n` and `n⋅∂ₙu = 0` at a
no-slip wall — but the discrete gradient only inherits that cancellation if the
discrete boundary actually resolves the curvature. It did not, and the resulting
FD-vs-adjoint gap was δ-independent and concentrated at the leading edge (95% of
the CD gradient error in L2 sat in the six parameters with x ≤ 0.15).

Note that this is a *global* option — Gmsh has no per-entity form of it — so it
also refines the farfield inlet arc (κ = 1/H), which is not wanted. The inlet is
therefore pinned transfinite in [`build_unstructured_msh`](@ref); see the comment
there.
"""
function set_curvature_sizing!(θ_max_deg::Real, h_min::Real)
    if θ_max_deg > 0
        gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", 360 / θ_max_deg)
        gmsh.option.setNumber("Mesh.MeshSizeMin", h_min)
    else
        gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", 0)
    end
    return nothing
end

"""
    build_unstructured_msh(am, airfoil_design, iter, chord, folder)

Build an unstructured mesh with a boundary layer around the airfoil. Must run
inside a Gmsh session (see [`with_gmsh`](@ref)); returns the written `.msh` file path.

Surface spacing is `rpoint` away from high-curvature regions and
curvature-limited near the nose — see [`set_curvature_sizing!`](@ref) and the
`curv_turn_deg` field of [`MeshSize`](@ref).
"""
function build_unstructured_msh(am::AirfoilMesh, airfoil_design::AirfoilDesign, iter::Int64, chord::Real, folder::String)

    function split_splines_points(airfoil_points::AirfoilPoints, AoA::Float64; pos=0.065, chord = 1.0)
         
        @unpack xu,xl,yu,yl = airfoil_points
        xur,yur = zeros(length(xu)),zeros(length(xu))
    
        for (i,(x,y)) in enumerate(zip(xu,yu))
            xur[i],yur[i] = rotate_points([x,y], AoA)
        end
    
        xlr,ylr = zeros(length(xl)),zeros(length(xl))
        for (i,(x,y)) in enumerate(zip(xl,yl))
            xlr[i],ylr[i] = rotate_points([x,y], AoA)
        end
    
        X = vcat(xur, xlr)
        Y =  vcat(yur, ylr)
    
         return X,Y
    end

    

    airfoil_points = airfoil_design.ap
    @unpack  Lback, H, meshref,BL_fl,BL_tt, airfoil_divisions, curv_turn_deg = am.MS
    @unpack AoA = am

    # `meshref` scales the mesh DENSITY: every isotropic size below is divided by
    # it, so meshref = 1 reproduces exactly the sizes this generator used before
    # the knob was wired up, and meshref = 2 halves them. `meshref = 0` is the
    # documented "very coarse, useful to debug" setting and is mapped to a fixed
    # large factor rather than to 1/0.
    #
    # BL_fl and BL_tt are deliberately NOT scaled: A1_convergence sweeps
    # (meshref, BL_fl) as an explicit pair, and scaling the first-layer height
    # here too would refine the boundary layer twice per level.
    hscale = meshref > 0 ? 1.0 / meshref : 20.0

    rpoint = 0.005 * hscale        # airfoil surface point spacing, away from the nose
    size_far = 0.01 * hscale       # BL field far size
    box_in = 0.03 * hscale         # wake box, refined
    box_out = 0.5 * hscale         # background, everywhere else

    gmsh.model.add("Model1")
    Lback = Lback*chord
    H= H*chord

    gmsh.model.geo.addPoint(Lback, -H, 0)
    gmsh.model.geo.addPoint(Lback, H, 0)
    
    gmsh.model.geo.addPoint(0.0, -H, 0)


    
    gmsh.model.geo.addPoint(0.0, H, 0)
        
    
    airfoil_gmsh_points = Int32[]



    X,Y = split_splines_points(airfoil_points, AoA)

    # Uniform `rpoint` everywhere except where the geometry actually turns (the
    # nose), which is handled by the curvature size source set below — these
    # per-point sizes are only honoured at the two spline endpoints anyway.
    for (xp,yp) in zip(X,Y)
            idx = gmsh.model.geo.addPoint(xp, yp, 0, rpoint)
            push!(airfoil_gmsh_points,idx)
    end
    
    
    trailing_coordinate =rotate_points([1.0,0.0],AoA)
    trailing = gmsh.model.geo.addPoint(trailing_coordinate[1],trailing_coordinate[2], 0, rpoint)




    origin_point = gmsh.model.geo.addPoint(0.0, 0.0, 0)

    limits_lines = zeros(Int32,2)
    inlet_lines = zeros(Int32,1)
    outlet_lines = zeros(Int32,1)
    


    #External Boundary Lines
    limits_lines[1] =  gmsh.model.geo.addLine(2,4)
    limits_lines[2] =  gmsh.model.geo.addLine(1,3)


    outlet_lines[1] = gmsh.model.geo.addLine(2,1)

    inlet_lines[1] = gmsh.model.geo.addCircleArc(4,origin_point,3)

    


    #Airfoil Splines
   
    airfoil_line = zeros(Int32,3)
    airfoil_line[1] = gmsh.model.geo.addSpline(airfoil_gmsh_points )
    airfoil_line[2] =  gmsh.model.geo.addLine(airfoil_gmsh_points[end],trailing)
    airfoil_line[3] =  gmsh.model.geo.addLine(trailing,airfoil_gmsh_points[1])

    # #Curve Loops
    outerLoop = gmsh.model.geo.addCurveLoop([
        -limits_lines[1],
        outlet_lines[1],
        limits_lines[2],
        -inlet_lines[1],
    ])
    
    airfoilLoop = gmsh.model.geo.addCurveLoop(airfoil_line)
    
    gmsh.model.geo.addPlaneSurface([outerLoop,airfoilLoop ])

    # Pin the inlet arc to the background spacing.
    #
    # `Mesh.MeshSizeFromCurvature` (set below, for the airfoil nose) is a GLOBAL
    # option: it refines every curved entity, and the inlet is a circle arc of
    # curvature 1/H. Its rule is "(360/θ) elements per 2π of turning", so the
    # half-arc gets (360/θ)/2 elements no matter how far away it is — 36 at
    # θ = 5°, 72 at θ = 2.5° — and since Gmsh takes the MIN over size sources,
    # `box_out` can never coarsen it back. Pushing the farfield out does not
    # help; shrinking it makes things worse, because κ = 1/H grows.
    #
    # A transfinite curve overrides every size source on that curve, so the
    # inlet is sized by `box_out` alone and `curv_turn_deg` only ever touches the
    # airfoil, which is what it is for. Only the inlet needs this: the outlet and
    # the top/bottom limits are straight, so curvature sizing never sees them,
    # and pinning the outlet would also cut off the wake-box refinement exactly
    # where the wake leaves the domain.
    inlet_nodes = max(2, round(Int, π * H / box_out) + 1)
    gmsh.model.geo.mesh.setTransfiniteCurve(inlet_lines[1], inlet_nodes, "Progression", 1.0)


    # gmsh.model.geo.mesh.setTransfiniteCurve(airfoil_line[1]  , 403, "Bump", 1.2)
    # gmsh.model.geo.mesh.setTransfiniteCurve(airfoil_line[2]  , 2, "Progression", 1.0)
    # gmsh.model.geo.mesh.setTransfiniteCurve(airfoil_line[3]  , 2, "Progression", 1.0)

    
    gmsh.model.mesh.field.add("BoundaryLayer", 1)
    gmsh.model.mesh.field.setNumbers(1, "CurvesList", [airfoil_line...])
    gmsh.model.mesh.field.setNumber(1, "Size", BL_fl)    # first layer height
    gmsh.model.mesh.field.setNumber(1, "SizeFar", size_far)
    gmsh.model.mesh.field.setNumber(1, "Thickness", BL_tt)    # total thickness
    gmsh.model.mesh.field.setNumber(1, "Ratio", 1.12)          # growth rate
    gmsh.model.mesh.field.setNumber(1, "Quads", 0) # unstructured meshes are TRI-only: a quad BL inside a tri mesh is a mixed mesh GridapGmsh cannot read
    
    # gmsh.model.mesh.field.setNumbers(1, "FanPointsList", [trailing])
    # gmsh.option.setNumber("Mesh.BoundaryLayerFanElements", 11)

    
    gmsh.model.mesh.field.setAsBoundaryLayer(1)


    # Wake-refinement Box field (Field[2]); the ONLY background size field.
    # The BoundaryLayer field must NOT be part of the background sizing: it is
    # applied through setAsBoundaryLayer above, and evaluating it as a size field
    # is Gmsh-version-dependent — on 4.14-git its SizeFar leaks over the whole
    # domain, inflating the mesh from ~36k to ~4.9M elements.
    gmsh.model.mesh.field.add("Box", 2)
    gmsh.model.mesh.field.setNumber(2, "VIn", box_in)
    gmsh.model.mesh.field.setNumber(2, "VOut", box_out)
    gmsh.model.mesh.field.setNumber(2, "XMin", 0.9)
    gmsh.model.mesh.field.setNumber(2, "XMax", 3.0)
    gmsh.model.mesh.field.setNumber(2, "YMin", -0.25)
    gmsh.model.mesh.field.setNumber(2, "YMax", 0.35)

    gmsh.model.mesh.field.setAsBackgroundMesh(2)

    # Curvature-limited sizing on the airfoil nose (see `set_curvature_sizing!`).
    # Gmsh takes the MIN over every size source, so this only ever refines.
    set_curvature_sizing!(curv_turn_deg, rpoint / 32)


    gmsh.model.geo.synchronize()


    #Points
    
    gmsh.model.addPhysicalGroup(0, [airfoil_gmsh_points[end],airfoil_gmsh_points[1],trailing], -1, "airfoil")
    gmsh.model.addPhysicalGroup(0, [1,2,3,4],-1,"limits")
    
    #Lines
    gmsh.model.addPhysicalGroup(1, airfoil_line,-1, "airfoil")
    gmsh.model.addPhysicalGroup(1, limits_lines,-1, "limits")
    gmsh.model.addPhysicalGroup(1, outlet_lines,-1, "outlet")
    gmsh.model.addPhysicalGroup(1, inlet_lines,-1, "inlet")
   

    # #Surfaces
    gmsh.model.addPhysicalGroup(2,[1],-1, "fluid")
    
    gmsh.model.geo.synchronize()

    return generate_and_write_msh(folder, iter)
end
