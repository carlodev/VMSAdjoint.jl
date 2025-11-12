using Gmsh
import Gmsh: gmsh



# --- Shared utilities ---------------------------------------------------------
function rotate_points(v::Vector{Float64}, AoA::Float64)
    x, y = v
    xrt = x * cosd(AoA) + y * sind(AoA)
    yrt = -x * sind(AoA) + y * cosd(AoA)
    return xrt, yrt
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


function find_origin_idx(leading_edge_points::Vector)
    _, idx = findmin(norm.(leading_edge_points))
    return idx
end


function create_msh(am::AirfoilMesh, airfoil_design::AirfoilDesign,  pp::PhysicalParameters ; iter::Int64= 0)
    @unpack folder = am
    return create_msh(am::AirfoilMesh, airfoil_design::AirfoilDesign,pp, iter, pp.c,am.folder)
end



function create_msh(am::AirfoilMesh, airfoil_design::AirfoilDesign,  pp::PhysicalParameters, folder::String; iter::Int64= 0 , )
    return create_msh(am::AirfoilMesh, airfoil_design::AirfoilDesign, pp, iter, pp.c, folder)
end




# --- STRUCTURED mesh workflow -------------------------------------------------

include(joinpath(@__DIR__,"GmshStructured.jl"))

function create_msh(am::AirfoilMesh{Structured},
                    airfoil_design::AirfoilDesign,
                    pp::PhysicalParameters, iter::Int, chord::Real, folder::String)
    return create_structured_msh(am, airfoil_design, iter, chord, folder)
end

# --- UNSTRUCTURED mesh workflow -----------------------------------------------

include(joinpath(@__DIR__,"GmshUnstructured.jl"))

function create_msh(am::AirfoilMesh{Unstructured},
                    airfoil_design::AirfoilDesign,
                    pp::PhysicalParameters, iter::Int, chord::Real, folder::String)
    return create_unstructured_msh(am, airfoil_design, iter, chord, folder)
end


