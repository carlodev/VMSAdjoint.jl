module Equations
using Gridap

using Parameters
using VMSAdjoint.Interfaces
using SegregatedVMSSolver.ParametersDef
using SegregatedVMSSolver.Equations

export equations_primal


export eq_direct_differentiation_steady
export eq_direct_differentiation_unsteady

include("PrimalEquations.jl")
# include("AdjointEquations.jl")
# include("DirectDifferentiation.jl")

end