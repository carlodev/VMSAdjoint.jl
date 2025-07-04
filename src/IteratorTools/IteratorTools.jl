module IteratorTools

using Parameters
using Gridap
using GridapGmsh
using Gridap.FESpaces
using Gridap.CellData:OperationCellField, GenericMeasure
using SegregatedVMSSolver.ParametersDef
using AirfoilTools
using VMSAdjoint.ParametersAdj
using VMSAdjoint.Interfaces
using VMSAdjoint.IncompressibleSolvers
using VMSAdjoint.IncompressibleSolvers: create_primal_spaces,create_adjoint_spaces,solve_inc_primal_steady
using JLD2
# using NLopt
using Optim
using LineSearches

using Interpolations
using ForwardDiff

export compute_airfoil_coefficients
export obj_fun
export dJobj_fun
export thickness_penalty
include("ObjectiveFunctions.jl")

export regularize_airfoil
include("RegularizationTools.jl")

export compute_sensitivity
include("Sensitivity.jl")

export solve_adjoint_optimization
include("OptLoop.jl")

export finite_difference_analysis
include("FDTools.jl")



end