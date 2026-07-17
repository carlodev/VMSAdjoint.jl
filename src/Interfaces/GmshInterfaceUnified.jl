using Gmsh
import Gmsh: gmsh

#############################################################################
# Shared geometric utilities
#############################################################################

"Rotate the point `(x, y)` by the angle of attack `AoA` (degrees, clockwise)."
function rotate_points(v::Vector{Float64}, AoA::Float64)
    x, y = v
    xrt = x * cosd(AoA) + y * sind(AoA)
    yrt = -x * sind(AoA) + y * cosd(AoA)
    return xrt, yrt
end

"Rotate the point set `[xs, ys]` by `AoA`, returning `[xr, yr]`."
function rotate_points(Mpoints::Vector{Vector{Float64}}, AoA::Float64)
    rotated = [rotate_points([x, y], AoA) for (x, y) in zip(Mpoints...)]
    return [first.(rotated), last.(rotated)]
end

function find_origin_idx(leading_edge_points::Vector)
    _, idx = findmin(norm.(leading_edge_points))
    return idx
end

#############################################################################
# Gmsh session helpers
#############################################################################

"""
    with_gmsh(f)

Run `f()` inside a fresh Gmsh session, guaranteeing `gmsh.finalize()` even when
`f` throws. A stale session left by a previously failed run is torn down first,
so mesh-generation failures never leak state into the next attempt.
"""
function with_gmsh(f::Function)
    gmsh.isInitialized() == 1 && gmsh.finalize()
    gmsh.initialize()
    try
        return f()
    finally
        gmsh.finalize()
    end
end

"""
    force_all_quads!()

Recombine the transfinite surfaces of a structured mesh into a pure-quad mesh
(GridapGmsh requires a single 2D cell type). Only meaningful for
`AirfoilMesh{Structured}` with `elements = :QUAD`; unstructured meshes are
always triangular (see `validate_mesh_config`).
"""
function force_all_quads!()
    gmsh.option.setNumber("Mesh.RecombineAll", 1)
end

"Generate the 2D mesh and write it to `<folder>/Mesh<iter>.msh`, returning the file path."
function generate_and_write_msh(folder::String, iter::Int)
    mkpath(folder)
    mesh_filename = joinpath(folder, "Mesh$iter.msh")
    gmsh.model.mesh.generate(2)
    gmsh.write(mesh_filename)
    return mesh_filename
end

#############################################################################
# Mesh-configuration validation
#############################################################################

allowed_elements(::AirfoilMesh{Structured})   = (:TRI, :QUAD)
allowed_elements(::AirfoilMesh{Unstructured}) = (:TRI,)

"""
    validate_mesh_config(am::AirfoilMesh, airfoil_design::AirfoilDesign)

Enforce the supported mesh/design combinations:
- `Structured` meshes support `:TRI` and `:QUAD` elements; `Unstructured` meshes
  only `:TRI` (a recombined unstructured mesh leaves mixed cell types that
  GridapGmsh cannot read).
- CST designs require a `Structured` mesh: the CST sensitivity kernel
  (`deformation_normal_field`) subtracts nodal position vectors between the
  baseline and perturbed models DOF-by-DOF, which is only well-defined when
  remeshing preserves the mesh topology — the transfinite structured mesh
  guarantees this, the unstructured mesher does not.
"""
function validate_mesh_config(am::AirfoilMesh{S}, airfoil_design::AirfoilDesign) where {S<:MeshStructure}
    if am.elements ∉ allowed_elements(am)
        throw(ArgumentError("elements = :$(am.elements) is not supported for AirfoilMesh{$S}; allowed: $(allowed_elements(am))"))
    end
    if airfoil_design isa AirfoilCSTDesign && !(am isa AirfoilMesh{Structured})
        throw(ArgumentError("CST designs require AirfoilMesh{Structured}: the perturbed geometry must be remeshed with identical topology for the DOF-wise mesh-velocity computation"))
    end
    return nothing
end

#############################################################################
# Entry points, dispatched on the mesh structure
#############################################################################

"""
    create_msh(am::AirfoilMesh, airfoil_design, pp::PhysicalParameters, [folder]; iter=0)

Build the Gmsh mesh for `airfoil_design` and return the `.msh` file path.
The mesh topology (structured C-type vs unstructured with boundary layer) is
selected by the `AirfoilMesh{Structured|Unstructured}` type parameter; the
element type by `am.elements` (see [`validate_mesh_config`](@ref) for the
allowed combinations).
"""
function create_msh(am::AirfoilMesh, airfoil_design::AirfoilDesign, pp::PhysicalParameters; iter::Int64=0)
    return create_msh(am, airfoil_design, pp, iter, pp.c, am.folder)
end

function create_msh(am::AirfoilMesh, airfoil_design::AirfoilDesign, pp::PhysicalParameters, folder::String; iter::Int64=0)
    return create_msh(am, airfoil_design, pp, iter, pp.c, folder)
end

include(joinpath(@__DIR__, "GmshStructured.jl"))

function create_msh(am::AirfoilMesh{Structured},
    airfoil_design::AirfoilDesign,
    pp::PhysicalParameters, iter::Int, chord::Real, folder::String)
    validate_mesh_config(am, airfoil_design)
    return with_gmsh() do
        build_structured_msh(am, airfoil_design, iter, chord, folder)
    end
end

include(joinpath(@__DIR__, "GmshUnstructured.jl"))

function create_msh(am::AirfoilMesh{Unstructured},
    airfoil_design::AirfoilDesign,
    pp::PhysicalParameters, iter::Int, chord::Real, folder::String)
    validate_mesh_config(am, airfoil_design)
    return with_gmsh() do
        build_unstructured_msh(am, airfoil_design, iter, chord, folder)
    end
end
