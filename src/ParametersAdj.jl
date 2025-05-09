using SegregatedVMSSolver
using Parameters
using SegregatedVMSSolver.ParametersDef
using AirfoilTools.Core
using AirfoilTools.AirfoilCST

export AirfoilMesh

struct AdjBC
    tagname::String #for multiple tags, we can put this as Vector{String}
    bc::Vector #Adjoint BC on the tagname; it depends on what we want to optimize
end



@with_kw struct AirfoilMesh <:MeshInfo
    AoA::Real #Angle of Attack - degrees
    cstdesign::AirfoilCSTDesign
    meshref::Int64=1.0
end