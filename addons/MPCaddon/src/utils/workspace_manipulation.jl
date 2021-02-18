# @Author: Massimo De Mauri <massimo>
# @Date:   2020-12-29T17:43:00+01:00
# @Email:  massimo.demauri@protonmail.com
# @Filename: mpc_workspace_manipulation.jl
# @Last modified by:   massimo
# @Last modified time: 2021-02-12T14:49:02+01:00
# @License: LGPL-3.0
# @Copyright: {{copyright}}


## Branch and Bound

# this function removes the constraits relative to the first part of the mpc short horizon problem
function mpc_remove_constraints!(workspace::BBworkspace,numSteps::Int,numVarsPerStep::Int;
                                 suppressErrors::Bool=false,localOnly::Bool=false)::Nothing

    @sync if !localOnly && !(workspace.sharedMemory isa NullSharedMemory)
         # call the local version of the function on the remote workers
         for p in workers()[1:workspace.settings.numProcesses-1]
             @async remotecall_fetch(OpenBB.eval,p,:(
                OpenBB.mpc_remove_constraints!(workspace,$numSteps,$numVarsPerStep;suppressErrors=$suppressErrors,localOnly=true)
             ))
         end

         # call the local version of the function on the main process
         mpc_remove_constraints!(workspace,numSteps,numVarsPerStep;suppressErrors=suppressErrors,localOnly=true)

    else

        # find to which time step each constraint refers
        timeSteps = div.(get_firstNZs(workspace.problem.cnsSet,2).+(numVarsPerStep-1),numVarsPerStep)
        if !issorted(timeSteps)
            permutation = sortperm(timeSteps,alg=MergeSort) # merge sort for stability (not the fastest)
            permute!(timeSteps,permutation)
            permute_constraints!(workspace,permutation,localOnly=true)
        end

        # remove the old constraints and the terminal ones
        indicesCnssToRemove = vcat(collect(1:findfirst((x)->(x==numSteps),timeSteps)-1), # eliminate first constraints
								   collect(findlast((x)->(x<timeSteps[end]),timeSteps):length(timeSteps)) # eliminate last constraints
								   )
        remove_constraints!(workspace,indicesCnssToRemove,suppressErrors=true,localOnly=true)

        # shift the remaining constraints backward
        remove_variables!(workspace.problem.cnsSet,collect(1:numVarsPerStep*numSteps))
        append_variables!(workspace.problem.cnsSet,numVarsPerStep*numSteps)

    end

    return
end



## Hybrid Branch and Bound
