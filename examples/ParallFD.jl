using ForwardDiff, MPI, PartitionedArrays

j(x) = x[1]^2 + x[2]^3


function dj(fun::Function, CDCL::Vector{Float64})
    @assert length(CDCL) == 2
    ForwardDiff.gradient(x -> fun(x), CDCL)
end


dj(j, [1.0,2.0])