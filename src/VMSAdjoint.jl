module VMSAdjoint


include(joinpath("ParametersAdj.jl","ParametersAdj.jl"))
include(joinpath("Interfaces","Interfaces.jl"))
include(joinpath("Equations","Equations.jl"))
include(joinpath("Iterators","Iterators.jl"))

include(joinpath("IncompressibleSolvers","IncompressibleSolvers.jl"))


# Re-export everything from submodules
include("Export.jl")


end
