module Interfaces

using LinearAlgebra
using Parameters
using GridapDistributed
using PartitionedArrays
using MPI
using SegregatedVMSSolver

using Gridap

using AirfoilTools
using Gmsh
using GridapGmsh
using VMSAdjoint.ParametersAdj
using SegregatedVMSSolver.ParametersDef
using SegregatedVMSSolver.Equations
using SegregatedVMSSolver.ExportUtility: unwrap_vector, wrap_vector, conv_VectorValue, get_dimension, conv_to_df,get_visgrid


export updatekey
export verifykey
include("ParamsInterfaces.jl")

export create_msh
include("GmshInterface.jl")

export get_aerodynamic_features
include("GridapInterface.jl")


include("VMSInterface.jl")


end