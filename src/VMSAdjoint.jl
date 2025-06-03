module VMSAdjoint

# Common packages used by all modules
using LinearAlgebra
using Parameters
using AirfoilTools
using Gridap
using GridapDistributed
using PartitionedArrays
using MPI
using SegregatedVMSSolver

include(joinpath("ParametersAdj","ParametersAdj.jl"))
include(joinpath("Interfaces","Interfaces.jl"))
include(joinpath("Equations","Equations.jl"))
include(joinpath("IncompressibleSolvers","IncompressibleSolvers.jl"))
include(joinpath("IteratorTools","IteratorTools.jl"))
include(joinpath("Main.jl"))



# Re-export everything from submodules
include("Export.jl")


end
