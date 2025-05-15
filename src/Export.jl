# Bring submodules into scope
using .ParametersAdj
using .Interfaces
using .Equations
using .IncompressibleSolvers
using .Iterators

# Automatically re-export what they export
for mod in (ParametersAdj, Interfaces, Equations, IncompressibleSolvers, Iterators)
    for name in names(mod, all = false, imported = false)
        @eval export $name
    end
end