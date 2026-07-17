"""
    Interfaces

Interfaces with the external packages: Gmsh (mesh generation), Gridap
(boundary/field extraction) and SegregatedVMSSolver (VMS stabilization).
"""
module Interfaces

using AirfoilTools
using Gmsh
using Parameters
using LinearAlgebra
using VMSAdjoint.ParametersAdj
using SegregatedVMSSolver.ParametersDef
using SegregatedVMSSolver.Equations

export updatekey
export verifykey
include("ParamsInterfaces.jl")

export create_msh
include("GmshInterfaceUnified.jl")

export get_aerodynamic_features
include("GridapInterface.jl")

include("VMSInterface.jl")

end
