
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
    gmsh.model.geo.addPoint(chord + slant , -H, 0)
    
    pback = rotate_points([Lback,0.0], AoA)
    gmsh.model.geo.addPoint(Lback, pback[2], 0)

    
    gmsh.model.geo.addPoint(chord + slant, H, 0)
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

    trailing_e_top = gmsh.model.geo.addPoint(trailing_coordinate[1] + slant,trailing_coordinate[2]+offset, 0)
    trailing_e_bottom = gmsh.model.geo.addPoint(trailing_coordinate[1] + slant,trailing_coordinate[2]-offset, 0)


    top_le_point = leading_edge_points[end]
    bottom_le_point = leading_edge_points[1]

    top_e_le_point = gmsh.model.geo.addPoint(0, offset, 0)
    bottom_e_le_point = gmsh.model.geo.addPoint(0,-offset, 0)

    top_back_point = gmsh.model.geo.addPoint(Lback, pback[2]+offset, 0)
    bottom_back_point = gmsh.model.geo.addPoint(Lback, pback[2]-offset, 0)

    origin_point0 = find_origin_idx(leading_edge_points_coordinates)
    origin_point = gmsh.model.geo.addPoint(0.0, 0.0, 0)

    limits_lines = zeros(Int32,4)
    inlet_lines = zeros(Int32,1)
    outlet_lines = zeros(Int32,4)
    


    #External Boundary Lines
    limits_lines[1] =  gmsh.model.geo.addLine(3,4)
    limits_lines[2] =  gmsh.model.geo.addLine(4,1)
    limits_lines[3] =  gmsh.model.geo.addLine(7,6)
    limits_lines[4] =  gmsh.model.geo.addLine(6,2)


    outlet_lines[1] = gmsh.model.geo.addLine(bottom_back_point,1)
    outlet_lines[2] = gmsh.model.geo.addLine(5,bottom_back_point)
    outlet_lines[3] = gmsh.model.geo.addLine(5,top_back_point)
    outlet_lines[4] = gmsh.model.geo.addLine(top_back_point,2)
    inlet_lines[1] = gmsh.model.geo.addCircleArc(7,origin_point,3)

    
    #Internal Lines
    gmsh.model.geo.addLine(trailing,5)
    gmsh.model.geo.addLine(trailing,trailing_e_top)
    gmsh.model.geo.addLine(trailing,trailing_e_bottom)

    gmsh.model.geo.addLine(trailing_e_top,6)
    gmsh.model.geo.addLine(trailing_e_top,top_back_point)

    gmsh.model.geo.addLine(trailing_e_bottom,4)
    gmsh.model.geo.addLine(trailing_e_bottom,bottom_back_point)

    gmsh.model.geo.addLine(top_e_le_point,trailing_e_top)
    gmsh.model.geo.addLine(bottom_e_le_point,trailing_e_bottom)

    gmsh.model.geo.addCircleArc(top_e_le_point,origin_point,bottom_e_le_point)

    gmsh.model.geo.addLine(top_e_le_point,7)
    gmsh.model.geo.addLine(bottom_e_le_point,3)

    gmsh.model.geo.addLine(top_le_point,top_e_le_point)
    gmsh.model.geo.addLine(bottom_le_point,bottom_e_le_point)

    #Airfoil Splines
   
    top_spline = gmsh.model.geo.addSpline([top_le_point, top_points...,trailing])
    bottom_spline = gmsh.model.geo.addSpline([bottom_le_point,bottom_points...,trailing])
    gmsh.model.geo.addSpline(leading_edge_points)
    
    

    #Curve Loops
    gmsh.model.geo.addCurveLoop([9,-21,-19,20])
    gmsh.model.geo.addCurveLoop([19,-23,26,22])
    gmsh.model.geo.addCurveLoop([21,1,-15,-18])
    gmsh.model.geo.addCurveLoop([-17,20,3,-13])
    gmsh.model.geo.addCurveLoop([-18,-23,25,12])
    gmsh.model.geo.addCurveLoop([-22,24,11,-17])
    gmsh.model.geo.addCurveLoop([15,2,-5,-16])
    gmsh.model.geo.addCurveLoop([12,16,-6,-10])
    gmsh.model.geo.addCurveLoop([10,7,-14,-11])
    gmsh.model.geo.addCurveLoop([13,4,-8,-14])
    
    
    for i = 1:1:10
        gmsh.model.geo.addPlaneSurface([i])
    end
    
    
    #vertical outer lines
    for i in [20,13,8,21,15,5]
        gmsh.model.geo.mesh.setTransfiniteCurve(i, 40, "Progression", 1.02)
    end

    if mesh_ref ==0
        corr = 0.1
    else
        corr = 1.0
    end

        
    #vertical inner lines
    for i in [7,11,22,23,12,6]
        gmsh.model.geo.mesh.setTransfiniteCurve(i, maximum([Int32(60*corr), Int32(round(70*mesh_ref))]), "Progression", 1.05) # 1.02
    end
    
    
    #inlet and leading edge
    for i in [9,19,26]
        gmsh.model.geo.mesh.setTransfiniteCurve(i, maximum([Int32(50*corr), Int32(round(40*mesh_ref))]), "Progression", 1.0)

    end

    #top airfoil
    for i in [24,17,3]
        gmsh.model.geo.mesh.setTransfiniteCurve(i, maximum([Int32(100*corr), Int32(round(50*mesh_ref))]), "Progression", 1.0)
    end
    
    #bottom airfoil
    for i in [25,18,1]
        gmsh.model.geo.mesh.setTransfiniteCurve(i, maximum([Int32(100*corr), Int32(round(50*mesh_ref))]), "Progression", 1.0)
    end
    
 
    
    #Shear Curves
    for i in [4,14,10,16,2]
        gmsh.model.geo.mesh.setTransfiniteCurve(i, maximum([Int32(40*corr), Int32(round(30*mesh_ref))]), "Progression", 1.05)
    end
    
    
    for i = 1:1:10
        gmsh.model.geo.mesh.setTransfiniteSurface(i)
    end
    
    
    
    gmsh.model.geo.synchronize()
    gmsh.option.setNumber("Mesh.RecombineAll", 1)
    
    #Points
    gmsh.model.addPhysicalGroup(0, [trailing,top_le_point,bottom_le_point], -1, "airfoil")
    gmsh.model.addPhysicalGroup(0, [trailing], -1, "trailing")


    gmsh.model.addPhysicalGroup(0, [top_e_le_point], -1, "top_e_le_point")

    gmsh.model.addPhysicalGroup(0, [5,top_back_point,bottom_back_point],-1, "outlet")
    gmsh.model.addPhysicalGroup(0, [3,4,1,7,6,2],-1,"limits")
    
    #Lines
    gmsh.model.addPhysicalGroup(1, [24,25,26],-1, "airfoil")
    gmsh.model.addPhysicalGroup(1, [1,2,3,4],-1, "limits")
    gmsh.model.addPhysicalGroup(1, [5,6,7,8],-1, "outlet")
    gmsh.model.addPhysicalGroup(1, [9],-1, "inlet")
   

    #Surfaces
    gmsh.model.addPhysicalGroup(2, collect(1:10),-1, "fluid")
    
    mkpath(folder)

    mesh_filename = joinpath(folder,"Mesh$iter.msh")

    gmsh.model.mesh.generate(2)
    gmsh.write(mesh_filename)
    gmsh.finalize()
    return mesh_filename
end
