


function create_unstructured_msh(am::AirfoilMesh, airfoil_design::AirfoilDesign, iter::Int64, chord::Real, folder::String)

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
    @unpack  Lback, H, meshref,BL_fl,BL_tt, airfoil_divisions = am.MS
    @unpack AoA = am

    gmsh.initialize()
    
    rpoint = 0.005
    gmsh.model.add("Model1")
    Lback = Lback*chord
    H= H*chord

    gmsh.model.geo.addPoint(Lback, -H, 0)
    gmsh.model.geo.addPoint(Lback, H, 0)
    
    gmsh.model.geo.addPoint(0.0, -H, 0)


    
    gmsh.model.geo.addPoint(0.0, H, 0)
        
    
    airfoil_gmsh_points = Int32[]



    X,Y = split_splines_points(airfoil_points, AoA)
    
    

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

    # gmsh.model.geo.mesh.setTransfiniteCurve( limits_lines[2], 20, "Progression", 1.0)
    # gmsh.model.geo.mesh.setTransfiniteCurve(-limits_lines[1], 20, "Progression", 1.0)

    # gmsh.model.geo.mesh.setTransfiniteCurve(outlet_lines[1], 20, "Progression", 1.0)
    # gmsh.model.geo.mesh.setTransfiniteCurve(-inlet_lines[1]  , 40, "Progression", 1.0)
  
      
    # gmsh.model.geo.mesh.setTransfiniteCurve(airfoil_line[1]  , 403, "Bump", 1.2)
    # gmsh.model.geo.mesh.setTransfiniteCurve(airfoil_line[2]  , 2, "Progression", 1.0)
    # gmsh.model.geo.mesh.setTransfiniteCurve(airfoil_line[3]  , 2, "Progression", 1.0)

    
    gmsh.model.mesh.field.add("BoundaryLayer", 1)
    gmsh.model.mesh.field.setNumbers(1, "CurvesList", [airfoil_line...])
    gmsh.model.mesh.field.setNumber(1, "Size", BL_fl)    # first layer height
    gmsh.model.mesh.field.setNumber(1, "SizeFar", 0.01) 
    gmsh.model.mesh.field.setNumber(1, "Thickness", BL_tt)    # total thickness
    gmsh.model.mesh.field.setNumber(1, "Ratio", 1.12)          # growth rate
    gmsh.model.mesh.field.setNumber(1, "Quads", 1)  
    
    # gmsh.model.mesh.field.setNumbers(1, "FanPointsList", [trailing])
    # gmsh.option.setNumber("Mesh.BoundaryLayerFanElements", 11)

    
    gmsh.model.mesh.field.setAsBoundaryLayer(1)


    # Create the Box field (Field[2])
    gmsh.model.mesh.field.add("Box", 2)
    gmsh.model.mesh.field.setNumber(2, "VIn", 0.03)
    gmsh.model.mesh.field.setNumber(2, "VOut", 0.5)
    gmsh.model.mesh.field.setNumber(2, "XMin", 0.9)
    gmsh.model.mesh.field.setNumber(2, "XMax", 3.0)
    gmsh.model.mesh.field.setNumber(2, "YMin", -0.25)
    gmsh.model.mesh.field.setNumber(2, "YMax", 0.35)
    


    # If you already had Field[1] for boundary layer, you can combine:
    gmsh.model.mesh.field.add("Min", 3)
    gmsh.model.mesh.field.setNumbers(3, "FieldsList", [1, 2])

    # Set background mesh field
    gmsh.model.mesh.field.setAsBackgroundMesh(3)


    gmsh.model.geo.synchronize()


    gmsh.option.setNumber("Mesh.RecombineAll", 1)

    gmsh.model.geo.synchronize()

    # gmsh.option.setNumber("Mesh.SubdivisionAlgorithm", 1)  # or 0, depending on surface shape

    
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
    
    mkpath(folder)

    gmsh.model.geo.synchronize()

    mesh_filename = joinpath(folder,"Mesh$iter.msh")

    
    gmsh.model.mesh.generate(2)
    
    gmsh.write(mesh_filename)
    gmsh.finalize()
    return mesh_filename
end
