# @Author: Massimo De Mauri <massimo>
# @Date:   2019-08-13T17:42:17+02:00
# @Email:  massimo.demauri@gmail.com
# @Filename: insert_node.jl
# @Last modified by:   massimo
# @Last modified time: 2019-10-23T18:21:12+02:00
# @License: LGPL-3.0
# @Copyright: {{copyright}}



# insert a node into a queue
function insert_node!(queue::Array{BBnode,1},node::BBnode,
                      settings::BBsettings,status::BBstatus)::Nothing

    if node.reliable || settings.unreliablesPriority == 0
    # normal queue insertion
        insertionPoint = 1
        for i in length(queue):-1:1
            if expansion_priority_rule(settings.expansionPriorityRule,node,queue[i],status)
                insertionPoint = i+1
                break
            end
        end
        insert!(queue,insertionPoint,node)

    elseif settings.unreliablesPriority == -1
        # put the new unreliable nodes at the bottom of the activeQueue to deal with them later
        insert!(workspace.activeQueue,0,node)

    elseif settings.unreliablesPriority == 1
        # put the new unreliable nodes at the top of the activeQueue to try to fastly get rid of them
        append!(workspace.activeQueue,node)
    else
        @error "wrong priority setting for unreliable nodes"
    end

    return
end


# insert a node into the BBtree
function insert_node!(workspace::BBworkspace,node::BBnode)::Nothing

    if isnan(node.objVal) # no info on the node
        # insert the node as first
        push!(workspace.activeQueue,node)

    elseif node.objVal - node.objGap > workspace.status.objUpB - workspace.settings.primalTolerance # the node is suboptimal

        # in interactive mode the suboptimal nodes are stored
        if workspace.settings.interactiveMode
            push!(workspace.unactivePool,node)
        end

    elseif node.objVal - node.objGap > workspace.settings.objectiveCutoff - workspace.settings.primalTolerance # the node objective is greater than the cutoff

        # declare the cutoff active
        workspace.status.cutoffActive = true

        # in interactive mode the suboptimal nodes are stored
        if workspace.settings.interactiveMode
            push!(workspace.unactivePool,node)
        end

    elseif node.avgAbsFrac == 0.0 # a new solution has been found

        if node.reliable || workspace.settings.acceptUnreliableSolutions # the solution is reliable or the algorithm is set to accept unreliable solutions

            # declare the cutoff inactive
            workspace.status.cutoffActive = false

            # insert new solution into the solutionPool
            push!(workspace.solutionPool,node)

            # update the number of solutions found
            workspace.status.numSolutions += 1
            # update the objective upper bound
            workspace.status.objUpB = node.objVal

            # update the global objective upper bound and the number of solutions found
            if !(workspace.sharedMemory isa NullSharedMemory)
                workspace.sharedMemory.stats[1] +=1
                if workspace.status.objUpB < workspace.sharedMemory.objectiveBounds[end]
                    workspace.sharedMemory.objectiveBounds[end] = workspace.status.objUpB
                end
            end

            # recompute optimality gaps
            if workspace.status.objUpB == Inf || workspace.status.objLoB == -Inf
                workspace.status.absoluteGap = workspace.status.relativeGap = Inf
            else
                workspace.status.absoluteGap = workspace.status.objUpB - workspace.status.objLoB
                workspace.status.relativeGap = workspace.status.absoluteGap/(1e-10 + abs(workspace.status.objUpB))
            end

        else # the solution is not reliable

            # the unreliability of the node cannot be solved... duly noted!
            workspace.status.reliable = false

            # store the unreliable solution in the unactivePool
            push!(workspace.unactivePool,node)
        end

    else # default behaviour
        insert_node!(workspace.activeQueue,node,workspace.settings,workspace.status)
    end

    return
end
