# @Author: Massimo De Mauri <massimo>
# @Date:   2019-05-20T10:04:17+02:00
# @Email:  massimo.demauri@gmail.com
# @Filename: update_settings.jl
# @Last modified by:   massimo
# @Last modified time: 2019-09-11T18:39:00+02:00
# @License: LGPL-3.0
# @Copyright: {{copyright}}


function update_objectiveCutoff!(workspace::BBworkspace{T1,T2,T3},newCutoff::Float64;
                                 suppressWarnings::Bool=false,localOnly::Bool=false)::Nothing where T1<:Problem where T2<:AbstractWorkspace where T3<:AbstractSharedMemory


    @sync if !localOnly && !(workspace.sharedMemory isa NullSharedMemory)
        # call the local version of the function on the remote workers
        for p in 2:workspace.settings.numProcesses
            @async remotecall_fetch(Main.eval,p,:(
                OpenBB.update_objectiveCutoff!(workspace,$newCutoff,suppressWarnings=$suppressWarnings,localOnly=true)
            ))
        end

        # call the function on the local worker
        update_objectiveCutoff!(workspace,newCutoff,suppressWarnings=suppressWarnings,localOnly=true)

    else
        # check the correctness of the input
        if !suppressWarnings && newCutoff > workspace.settings.objectiveCutoff && workspace.status.description != "new" && myid() == 1
            @warn "Relaxing the cutoff after some iterations may lead to incorrect results"
        end
        # change the cutoff
        workspace.settings.objectiveCutoff = newCutoff
        # invalitate the solutions that do not respect the new cutoff
        if workspace.status.objUpB > newCutoff
            # update the status
            workspace.status.objUpB = workspace.status.absoluteGap = workspace.status.relativeGap = Inf
            workspace.status.numSolutions = 0
            if workspace.status.description == "optimalSolutionFound"
                workspace.status.description = "interrupted"
            end

            append!(workspace.unactivePool,workspace.solutionPool)
            deleteat!(workspace.solutionPool,1:length(workspace.solutionPool))
        else
            solutionsToElim = Array{Int,1}()
            for k in 1:length(workspace.solutionPool)
                if workspace.solutionPool[k].objUpB > newCutoff - workspace.settings.primalTolerance
                    push!(solutionsToElim,k)
                    push!(workspace.unactivePool,workspace.solutionPool[k])
                    workspace.status.numSolutions -= 1
                    if !(workspace.sharedMemory isa NullSharedMemory)
                        workspace.sharedMemory.stats[1] -= 1
                    end
                end
            end
            deleteat!(workspace.solutionPool,solutionsToElim)
        end
    end

    # mark the workspace as outdated
    make_outdated!(workspace)

    return
end




#TODO: implement update_settings!
# function update_settings!(workspace::BBworkspace{T1,T2,T3};bb_settings::AbstractSettings=NullSettings(),ss_settings::AbstractSettings=NullSettings())::Nothing
#          where T1<:Problem where T2<:AbstractWorkspace where T3<:AbstractSharedMemory
#
#
#     @error "Not implemented yet"
#
#     # check if it is possible to make changes
#     @assert workspace.settings.interactiveMode == true || workspace.status.description == "new"
#
#
#     # update the BB settings
#     if !(bb_settings isa NullSettings)
#         if workspace.status.totalTime > 0.
#
#             # adapt the workspace to the new settings
#
#
#             # pass from interactive mode to static mode
#             if workspace.interactiveMode == true && bb_settings.interactiveMode == false
#                 for k in length(workspace.solutionPool):-1:1
#                     if workspace.solutionPool[k].reliable
#                         deleteat!(workspace.solutionPool,k+1:length(workspace.solutionPool))
#                         deleteat!(workspace.solutionPool,1:k-1)
#                         break
#                     elseif k == 1
#                         deleteat!(workspace.solutionPool,1:length(workspace.solutionPool)-1)
#                     end
#                 end
#                 deleteat!(workspace.unactivePool,1:length(workspace.unactivePool))
#             end
#         end
#
#         for field in fieldnames(BBsettings)
#             setfield!(workspace.settings,field,getfield(bb_settings,field))
#         end
#     end
#
#     # update the subsolver settings
#     if !(ss_settings isa NullSettings)
#         update_settings!(workspace.subsolverWS,ss_settings,bb_primalTolerance=workspace.settings.primalTolerance,bb_timeLimit=workspace.settings.timeLimit)
#     end
#
#     return
# end
