module VMSAdjoint


include(joinpath("ParametersAdj","ParametersAdj.jl"))
include(joinpath("Interfaces","Interfaces.jl"))
include(joinpath("Equations","Equations.jl"))
include(joinpath("IncompressibleSolvers","IncompressibleSolvers.jl"))

include(joinpath("IteratorTools","IteratorTools.jl"))



# Re-export everything from submodules
include("Export.jl")


end
