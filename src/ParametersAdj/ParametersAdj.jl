module ParametersAdj

"""
    In this module are defined the main struct used in the package. It is extending the basic ones from SegregatedVMSSolver.ParametersDef and AirfoilTools
"""

using SegregatedVMSSolver
using Parameters
using Gridap
using SegregatedVMSSolver.ParametersDef
using AirfoilTools

export AirfoilNormals
export AirfoilScalar
export AirfoilModel
export AirfoilMesh
export DesignBounds
export ThickPenalty
export AdjSolver
export AdjointProblem
export AdjSolution

include("ParamsCore.jl")


end