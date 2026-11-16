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
using Interpolations
using Gridap
using Gridap.Algebra
using Gridap.FESpaces
using Gridap.Geometry
using Gridap.CellData: get_cell_measure
using VMSAdjoint.ParametersAdj
using SegregatedVMSSolver.ParametersDef
using SegregatedVMSSolver.Equations

export updatekey
export verifykey
include("ParamsInterfaces.jl")

export create_msh
export rotate_points
include("GmshInterfaceUnified.jl")

export get_aerodynamic_features
export airfoil_boundary!
include("GridapInterface.jl")

export deform_airfoil_mesh
export deform_model
export mesh_motion_cache
export mesh_motion_cache!
export solve_mesh_displacement
export airfoil_surface_displacement
include("MeshMotion.jl")

include("VMSInterface.jl")

end
