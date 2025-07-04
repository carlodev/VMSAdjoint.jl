module Equations
using Gridap

using Parameters
using VMSAdjoint.Interfaces
using SegregatedVMSSolver
using SegregatedVMSSolver.ParametersDef
using SegregatedVMSSolver.Equations

export equations_primal
export equations_adjoint


include("PrimalEquations.jl")
include("AdjointEquations.jl")

end