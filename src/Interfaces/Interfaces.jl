module Interfaces

using AirfoilTools
using Gmsh
using Parameters
using LinearAlgebra
using VMSAdjoint.ParametersAdj
using SegregatedVMSSolver.ParametersDef
using SegregatedVMSSolver.Equations
using AirfoilTools

export updatekey
export verifykey
include("ParamsInterfaces.jl")

export create_msh
include("GmshInterface.jl")

export get_aerodynamic_features
include("GridapInterface.jl")


export get_aerodynamic_features
include("VMSInterface.jl")

export perturb_DesignParameter
include("CSTInterface.jl")

end