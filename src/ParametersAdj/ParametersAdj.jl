module ParametersAdj


using Parameters
using AirfoilTools
using SegregatedVMSSolver.ParametersDef

export AirfoilNormals
export AirfoilScalar
export AirfoilModel
export AirfoilMesh
export AdjSolver
export AdjointProblem
export AdjSolution

include("ParamsCore.jl")


end