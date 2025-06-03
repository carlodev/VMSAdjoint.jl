using SegregatedVMSSolver
using SegregatedVMSSolver.ParametersDef
using PartitionedArrays
using MPI
using VMSAdjoint.ParametersAdj

function solve(adjp::AdjointProblem, backend)

    backend() do distribute
        if backend == with_mpi
            comm = MPI.COMM_WORLD
            #To avoid multiple printing of the same line in parallel
            if MPI.Comm_rank(comm) != 0
              redirect_stderr(devnull)
              redirect_stdout(devnull)
            end
           
        end

        @sunpack rank_partition,order = adjp.vbcase
        
        parts  = distribute(LinearIndices((prod(rank_partition),)))
        @info "parts created"
        solve_adjoint_optimization(adjp, parts)

    end

return true

end
