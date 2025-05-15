module ParametersAdj

using SegregatedVMSSolver
using Parameters
using Gridap
using SegregatedVMSSolver.ParametersDef
using AirfoilTools

export AirfoilNormals
export AirfoilScalar
export AirfoilModel
export AdjBC
export AirfoilMesh
export AdjointProblem
export AdjSolution

include("ParamsCore.jl")


end