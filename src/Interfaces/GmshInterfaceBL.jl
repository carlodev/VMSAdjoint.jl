
using Gmsh
import Gmsh: gmsh

function split_splines_points(airfoil_points::AirfoilPoints, AoA::Float64; pos=0.065, chord = 1.0)
         
    @unpack xu,xl,yu,yl = airfoil_points

    threshold = pos * chord
    top_LE_point = findmin(abs.(xu .- threshold))[2]
    bottom_LE_point = findmin(abs.(xl .- threshold))[2]
   

    Mtop = rotate_points([xu[1:top_LE_point],yu[1:top_LE_point]], AoA)

    Mbottom = rotate_points([xl[bottom_LE_point:end],yl[bottom_LE_point:end]], AoA)
    
    Mle = rotate_points([ [xu[top_LE_point+1:end];xl[1:bottom_LE_point-1]],
    [yu[top_LE_point+1:end];yl[1:bottom_LE_point-1]]  ], AoA)

     return  reverse.(Mtop), Mbottom, reverse.(Mle)
end


function rotate_points(Mpoints::Vector{Vector{Float64}}, AoA::Float64)
    xr = Float64[]
    yr = Float64[]
    for (x,y) in zip(Mpoints...)
        xrt= x * cosd(AoA) + y*sind(AoA)
        yrt = -1*x * sind(AoA) + y *cosd(AoA)
        push!(xr,xrt)
        push!(yr,yrt)
    end
    return [xr,yr]
end

function rotate_points(v::Vector{Float64}, AoA::Float64)
    x = v[1] 
    y = v[2]
    
    xrt= x * cosd(AoA) + y*sind(AoA)
    yrt = -1*x * sind(AoA) + y *cosd(AoA)

    return xrt,yrt
end

function find_origin_idx(leading_edge_points::Vector)
    _,idx = findmin(norm.(leading_edge_points))
    return idx
end

function create_msh(am::AirfoilMesh, airfoil_design::AirfoilDesign,  pp::PhysicalParameters ; iter::Int64= 0)
    @unpack AoA, meshref, folder, H, Lback = am
    @unpack ap = airfoil_design
    return create_msh(ap;H=H, Lback =Lback, AoA=AoA, iter=iter, chord = pp.c, mesh_ref = meshref, folder = folder)
end



function create_msh(am::AirfoilMesh, airfoil_design::AirfoilDesign,  pp::PhysicalParameters, folder::String; iter::Int64= 0 , )
    @unpack AoA, meshref,H, Lback = am
    @unpack ap = airfoil_design
    return create_msh(ap; H=H, Lback =Lback, AoA=AoA, iter=iter, chord = pp.c, mesh_ref = meshref, folder = folder)
end

"""
    create_msh(airfoil_points::AirfoilPoints; AoA=0.0, iter = 0, chord= 1.0, mesh_ref=1.0)

From a set of `airfoil_points` it creates the .msh file. Incresing `mesh_ref` is increasing the mesh density.
"""
function create_msh(airfoil_points::AirfoilPoints; H=8.0, Lback =8.0, AoA=0.0, iter = 0, chord= 1.0, mesh_ref=1.0, folder="MeshFiles")


    gmsh.initialize()
    
    gmsh.model.add("Model1")
    Lback = Lback*chord
    H= H*chord
    offset = 2.35
    slant = 2.0
    
    gmsh.model.geo.addPoint(Lback, -H, 0)
    gmsh.model.geo.addPoint(Lback, H, 0)
    
    gmsh.model.geo.addPoint(0.0, -H, 0)
    # gmsh.model.geo.addPoint(chord + slant , -H, 0)
    
    # pback = rotate_points([Lback,0.0], AoA)
    # gmsh.model.geo.addPoint(Lback, pback[2], 0)

    
    # gmsh.model.geo.addPoint(chord + slant, H, 0)
    gmsh.model.geo.addPoint(0.0, H, 0)
    
    # gmsh.model.geo.addPoint(Lfront, 0.0, 0)
    
    
    top_points = Int32[]
    bottom_points = Int32[]
    leading_edge_points = Int32[]
    leading_edge_points_coordinates = Vector[]


    Mtop, Mbottom, Mle = split_splines_points(airfoil_points, AoA)
    
    

    for (xp,yp) in zip(Mtop...)
            idx = gmsh.model.geo.addPoint(xp, yp, 0)
            push!(top_points,idx)
    end
    
    for  (xp,yp) in zip(Mbottom...)
            idx = gmsh.model.geo.addPoint(xp, yp, 0)
            push!(bottom_points,idx)
    end
    
    for  (xp,yp) in zip(Mle...)
        idx = gmsh.model.geo.addPoint(xp, yp, 0)
        push!(leading_edge_points,idx)
        push!(leading_edge_points_coordinates, [xp,yp])
    end
    
    trailing_coordinate =rotate_points([1.0,0.0],AoA)

    trailing = gmsh.model.geo.addPoint(trailing_coordinate[1],trailing_coordinate[2], 0)


    top_le_point = leading_edge_points[end]
    bottom_le_point = leading_edge_points[1]



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
   
    top_spline = gmsh.model.geo.addSpline([top_le_point, top_points...,trailing])
    bottom_spline = gmsh.model.geo.addSpline([bottom_le_point,bottom_points...,trailing])
    leading_edge_spline= gmsh.model.geo.addSpline(leading_edge_points)
    
    

    # #Curve Loops
    gmsh.model.geo.addCurveLoop([- limits_lines[1],outlet_lines[1],  limits_lines[2], -inlet_lines[1]   ])

    gmsh.model.geo.addCurveLoop([- top_spline, -leading_edge_spline,  bottom_spline   ])


    gmsh.model.geo.addPlaneSurface([1,2])
    
    gmsh.model.geo.mesh.setTransfiniteCurve( limits_lines[2], 20, "Progression", 1.0)
    gmsh.model.geo.mesh.setTransfiniteCurve(-limits_lines[1], 20, "Progression", 1.0)

    gmsh.model.geo.mesh.setTransfiniteCurve(outlet_lines[1], 20, "Progression", 1.0)
    gmsh.model.geo.mesh.setTransfiniteCurve(-inlet_lines[1]  , 40, "Progression", 1.0)
  

    gmsh.model.geo.mesh.setTransfiniteCurve(top_spline  , 201, "Progression", 1.0)
    gmsh.model.geo.mesh.setTransfiniteCurve(bottom_spline  , 100, "Progression", 1.0)
    gmsh.model.geo.mesh.setTransfiniteCurve(leading_edge_spline  , 40, "Progression", 1.0)
    
    gmsh.model.mesh.field.add("BoundaryLayer", 1)
    gmsh.model.mesh.field.setNumbers(1, "CurvesList", [top_spline, bottom_spline,leading_edge_spline])
    gmsh.model.mesh.field.setNumber(1, "Size", 1e-4)    # first layer height
    gmsh.model.mesh.field.setNumber(1, "SizeFar", 0.01) 
    gmsh.model.mesh.field.setNumber(1, "Thickness", 0.02)    # total thickness
    gmsh.model.mesh.field.setNumber(1, "Ratio", 1.12)          # growth rate
    gmsh.model.mesh.field.setNumber(1, "Quads", 1)         
    gmsh.model.mesh.field.setAsBoundaryLayer(1)



    # Create the Box field (Field[2])
    gmsh.model.mesh.field.add("Box", 2)
    gmsh.model.mesh.field.setNumber(2, "VIn", 0.03)
    gmsh.model.mesh.field.setNumber(2, "VOut", 0.5)
    gmsh.model.mesh.field.setNumber(2, "XMin", 0.9)
    gmsh.model.mesh.field.setNumber(2, "XMax", 3.0)
    gmsh.model.mesh.field.setNumber(2, "YMin", -0.45)
    gmsh.model.mesh.field.setNumber(2, "YMax", 0.45)

    # If you already had Field[1] for boundary layer, you can combine:
    gmsh.model.mesh.field.add("Min", 3)
    gmsh.model.mesh.field.setNumbers(3, "FieldsList", [1, 2])

    # Set background mesh field
    gmsh.model.mesh.field.setAsBackgroundMesh(3)


    gmsh.model.geo.synchronize()
    gmsh.option.setNumber("Mesh.RecombineAll", 1)

    gmsh.option.setNumber("Mesh.SubdivisionAlgorithm", 1)  # or 0, depending on surface shape

    
    #Points
    # gmsh.model.addPhysicalGroup(0, [trailing,top_le_point,bottom_le_point, top_points..., bottom_points..., leading_edge_points...], -1, "airfoil")
    
    gmsh.model.addPhysicalGroup(0, [trailing,top_le_point,bottom_le_point,], -1, "airfoil")

    gmsh.model.addPhysicalGroup(0, [trailing], -1, "trailing")
    


    gmsh.model.addPhysicalGroup(0, [1,2,3,4],-1,"limits")
    
    #Lines
    gmsh.model.addPhysicalGroup(1, [top_spline, leading_edge_spline,  bottom_spline],-1, "airfoil")
    gmsh.model.addPhysicalGroup(1, [limits_lines...],-1, "limits")
    gmsh.model.addPhysicalGroup(1, [outlet_lines...],-1, "outlet")
    gmsh.model.addPhysicalGroup(1, [inlet_lines...],-1, "inlet")
   

    # #Surfaces
    gmsh.model.addPhysicalGroup(2,[1],-1, "fluid")
    
    mkpath(folder)

    mesh_filename = joinpath(folder,"Mesh$iter.msh")

    gmsh.model.mesh.generate(2)
    gmsh.write(mesh_filename)
    gmsh.finalize()
    return mesh_filename
end
