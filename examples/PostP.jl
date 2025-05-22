using Revise
using AirfoilTools
using VMSAdjoint
using Gridap, GridapGmsh
using SegregatedVMSSolver.ParametersDef
using Plots
using JLD2





adj_sol0 = load("results/ADJ_SOL0.jld2")["adj_sol"]
plot(adj_sol0.fgrad)