using Revise
using AirfoilTools
using AirfoilTools.AirfoilCST
using VMSAdjoint.Interfaces
using Gridap, GridapGmsh
AoA = 1.0
NW0 = 15 #number of w0 to parametrize top and bottom
cst0 = CST_NACA0012(;N=20, t=0.0)
cstd0 = AirfoilCSTDesign(cst0,200)


modelname =create_msh(cstd0.ap; AoA=AoA, mesh_ref=2.0)



model = GmshDiscreteModel(modelname)
writevtk(model, "model_original")

