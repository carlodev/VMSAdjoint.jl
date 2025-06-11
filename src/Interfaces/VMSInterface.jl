import SegregatedVMSSolver.Equations: G_params, compute_G, compute_GG, compute_gg

"""
  SegregatedVMSSolver.Equations.G_params

Extension of the function on single-processor
"""
function SegregatedVMSSolver.Equations.G_params(Ω::Triangulation, D::Int64)
    G = compute_G(Ω,D)
    GG = compute_GG(Ω,D)
    gg = compute_gg(Ω,D)
   return G, GG, gg
 end
 