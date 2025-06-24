# VMSAdjoint

[![Build Status](https://github.com/carlodev/VMSAdjoint.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/carlodev/VMSAdjoint.jl/actions/workflows/CI.yml?query=branch%3Amain)


## How to use the package:

## Tutorial
In this tutorial we are going to see how to obtain the sensitivities on a NACA 0012 airfoil using the Variational MultiScale Adjoint.

Start loading the packages
```julia
using AirfoilTools
using VMSAdjoint
using Gridap, GridapGmsh
using SegregatedVMSSolver.ParametersDef
```

Define the physical parameters of the simulation
```julia
fname = "n0012.csv" #file name with the coordinates to load
AoA = 2.5 #angle of attack, in degrees
Re = 1000 #Reynolds number
uin = [1.0,0.0] #inlet velocity
order = 2 #order of the elements
```

Define the parametrization of the airfoil. In this case we are using control points on the top and bottom surface of the airfoil. We define the x-coordinates of them, the y-coordinates are computed when we create the Design.

```julia
ap0 = get_airfoil_coordinates(joinpath(@__DIR__, fname))

#40 control points in total, 20 on the suction side, 20 on the pressure side
control_px = collect(LinRange(0.05,0.95,20))
control_points = ControlPoints(control_px,control_px)

R = 0.15 #support radius
RBFfun = RBFFunctionLocalSupport(RBF_CP4, R) #RBF_CP4 is the Radial Basis Function

rbfg = RBFGeometry(control_points,RBFfun)
rbfd = RBFDesign(rbfg, ap0)
```

We econde the values of the simulation.
```julia

sprob = StabilizedProblem(VMS(order))

physicalp = PhysicalParameters(Re=Re, u_in=uin)
timep = TimeParameters(dt=0.05, tF=1.0, time_window=(0.8, 1.0)) #the time-window define the time-span for time-averaging


meshinfo = AirfoilMesh(AoA= AoA, meshref=1) #increasing meshref, we increse the resolution of the mesh in the domain
meshp = MeshParameters((1,1), 2, meshinfo)
exportp = ExportParameters(printinitial=true,printmodel=true,name_tags=["airfoil"], fieldexport=[["uh","ph","friction"]])

solverp = SolverParameters(θ=1.0,  M=1, matrix_freq_update=4) 

simparams = SimulationParameters(timep,physicalp,solverp,exportp)

airfoil_case = Airfoil(meshp,simparams,sprob)
```

Now, the objective function is defined.  it only takes one argument [CD,CL], then you can define all the keywords you want. The boundary conditions are defined as -dJ/dCDCL

```julia
function J(CDCL; CLtarget=0.75)
    CD,CL=CDCL
    return CD/CL  
end  
```

Now we define the adjoint solver, where there are many customizable options, see in `src/ParametersAdj/ParametersAdj.jl`
```julia
adj_solver = AdjSolver(δ=0.0001)

```

Then, we solve the problem. `(:unsteady, :steady)` means that the primal flow is solved as `unsteady` and the adjoint is solve as `steady`.

```julia
adjoint_airfoil_problem = AdjointProblem( rbfd,airfoil_case,adj_solver, (:unsteady, :steady), J )

```




### Adjoint Simulation
Look at the example: EX_M.jl

### Finite Difference
Look at the example: EX_FD.jl





## Installation
The package has not been registered yet so you can install it as:
```julia
using Pkg
Pkg.add(url="https://github.com/carlodev/VMSAdjoint.jl")
```
### Packages recommendation

https://github.com/carlodev/SegregatedVMSSolver.jl#Update-deps
https://github.com/carlodev/AirfoilTools.jl#master
 
They will published soon - they are still under development

## Note
The package has been developed having in mind the parallelization, which is going to be implemented.

## Contributing
It is a collaborative project open to contributions. You can:
- Open a new issue
- Contact the project administator
- Open a PR with the contribution