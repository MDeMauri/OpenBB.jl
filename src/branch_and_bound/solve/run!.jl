# @Author: Massimo De Mauri <massimo>
# @Date:   2019-02-06T18:33:16+01:00
# @Email:  massimo.demauri@gmail.com
# @Project: OpenBB
# @Filename: run!.jl
# @Last modified by:   massimo
# @Last modified time: 2019-08-27T20:19:45+02:00
# @License: LGPL-3.0
# @Copyright: {{copyright}}


# branch and bound algorithm
function run!(workspace::BBworkspace{T1,T2})::Nothing where T1<:AbstractWorkspace where T2<:AbstractSharedMemory

    # timing
    lastTimeCheckpoint = time()

    # update algorithm status
    workspace.status.description = "running"

    # set the algorithm in "active state"
    idle = false

    # id of the process the workspace is being manipulated from
    processId = myid()

    # print the initial algorithm status
    printCountdown = 0.0

    # keep in memory how many nodes the algorithm has explored since the last time
    # it has sent a node to the neighbouring process
    timeToShareNodes = false

    # communicate the initial objective lower-bound to the other nodes
    if !(workspace.sharedMemory isa NullSharedMemory)
        # communicate the current local lower bound
        workspace.sharedMemory.objectiveBounds[processId] = workspace.status.objLoB
    end


    # main loop
    while !idle

        # get elapsed time for the last iteration
        iterationTime = time() - lastTimeCheckpoint
        lastTimeCheckpoint = time()
        # update total time
        workspace.status.totalTime += iterationTime
        printCountdown -= iterationTime

        # print algorithm status
        if workspace.settings.verbose && processId == 1
           if printCountdown <= 0
               print_status(workspace)
               printCountdown = workspace.settings.statusInfoPeriod
          end
        end

        # stopping conditions
        if length(workspace.activeQueue) == 0 || # no more nodes in the queue
           workspace.status.objLoB > workspace.settings.objectiveCutoff || # no solutions within the cutoff possible
           workspace.status.totalTime >= workspace.settings.timeLimit || # time is up
           workspace.settings.customStoppingRule(workspace) || # custom stopping rule triggered
           workspace.status.absoluteGap < workspace.settings.absoluteGapTolerance || # reached required absolute gap
           workspace.status.relativeGap < workspace.settings.relativeGapTolerance || # reached required relative gap
           (workspace.settings.numSolutionsLimit > 0 && workspace.status.numSolutions >= workspace.settings.numSolutionsLimit) # the required number of solutions has been found

           # set the algorithm in idle state
           if workspace.sharedMemory isa NullSharedMemory || !isready(workspace.sharedMemory.outputChannel)
               idle = true
           end

       else # continue with branch and bound

           # apply rounding heuristics
           if workspace.activeQueue[end].avgAbsFrac > 0 &&
              workspace.activeQueue[end].avgAbsFrac < workspace.settings.roundingHeuristicsThreshold &&
              (workspace.activeQueue[end].objective-workspace.status.objLoB)/(1e-10 + workspace.activeQueue[end].objective) < workspace.settings.roundingHeuristicsThreshold

               node = simple_rounding_heuristics(workspace.activeQueue[end],workspace)
               solve!(node,workspace)

               # new solution?
               if node.reliable &&
                  node.objective < min(workspace.status.objUpB,workspace.settings.objectiveCutoff) - workspace.settings.primalTolerance

                   # insert new solution into the solutionPool
                   push!(workspace.solutionPool,node)

                   # update the number of solutions found
                   workspace.status.numSolutions += 1
                   # update the objective upper bound
                   workspace.status.objUpB = node.objective

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
               end
           end

           # pick a node to process from the activeQueue
           node = pop!(workspace.activeQueue)


            # check if node is already suboptimal (the upper-bound might have changed)
            if node.objective > workspace.status.objUpB - workspace.settings.primalTolerance

                # in interactive mode the suboptimal nodes are stored
                if workspace.settings.interactiveMode
                    push!(workspace.unactivePool,node)
                end

            elseif node.objective > workspace.settings.objectiveCutoff - workspace.settings.primalTolerance

                # declare the cutoff active
                workspace.status.cutoffActive = true

                # in interactive mode the suboptimal nodes are stored
                if workspace.settings.interactiveMode
                    push!(workspace.unactivePool,node)
                end

            elseif !(workspace.sharedMemory isa NullSharedMemory) && timeToShareNodes && # we are multiprocessing and it is time to send a child to the next process
                   !isready(workspace.sharedMemory.outputChannel) # the channel is free

                # send the new node to the neighbouring process
                put!(workspace.sharedMemory.outputChannel,node)

                # the next node has to be explored locally
                timeToShareNodes = false

            else # explore the node locally

                # the next node will be sent to the neighbouring process
                timeToShareNodes = true

                # create a list of children nodes
                children = branch_and_solve!(node,workspace)

                # remove the infeasible children
                filter!((child)->child.objective<Inf,children)

                # insert the rest of the children in the BBtree (checking for solutions and suboptimals)
                for child in children
                    insert_node!(workspace,child)
                end
            end

            # recompute the objective lower bound
            if workspace.status.objLoB == -Inf || node.objective == workspace.status.objLoB || length(workspace.activeQueue) == 0
                # compute the new lower-bound
                newObjLoB = workspace.status.objUpB
                for i in length(workspace.activeQueue):-1:1
                    # if we have unreliable problems in the activeQueue we cannot update the lower bound
                    if !workspace.activeQueue[i].reliable
                        newObjLoB = workspace.status.objLoB
                        break
                    elseif workspace.activeQueue[i].objective < newObjLoB
                        newObjLoB = workspace.activeQueue[i].objective
                    end
                end

                if workspace.status.objLoB > newObjLoB + workspace.settings.primalTolerance
                    @warn "branch and bound: the objective lower bound has decreased from "*string(workspace.status.objLoB)*" to "*string(newObjLoB)*"..."
                end

                # update the objective lower-bound
                workspace.status.objLoB = newObjLoB

                # recompute optimality gaps
                if workspace.status.objUpB == Inf || workspace.status.objLoB == -Inf
                    workspace.status.absoluteGap = workspace.status.relativeGap = Inf
                else
                    workspace.status.absoluteGap = workspace.status.objUpB - workspace.status.objLoB
                    workspace.status.relativeGap = workspace.status.absoluteGap/(1e-10 + abs(workspace.status.objUpB))
                end
            end
        end

        # communicate with the other processes
        if !(workspace.sharedMemory isa NullSharedMemory)

            if idle
                # declare the worker ready to stop
                workspace.sharedMemory.arrestable[processId] = true

                # keep track of the time spent in waiting
                waitingStartTime = time()

                # wait for a message to arrive in the inputChannel
                while !isready(workspace.sharedMemory.inputChannel) &&
                      !all(@. workspace.sharedMemory.arrestable)
                end

                # keep track of the time spent in waiting
                workspace.status.waitingTime += time()-waitingStartTime
            end

            if isready(workspace.sharedMemory.inputChannel)
                if idle
                    # undeclare the worker ready to stop
                    workspace.sharedMemory.arrestable[processId] = false
                    # go active
                    idle = false
                end

                # take a new node from the input channel
                newNode = take!(workspace.sharedMemory.inputChannel)

                # insert the new node in the active queue
                insert_node!(workspace.activeQueue,node,workspace.settings,workspace.status)

                # update the objective lower bound
                if newNode.objective < workspace.status.objLoB
                    workspace.status.objLoB = newNode.objective
                    # recompute optimality gaps
                    if workspace.status.objUpB == Inf || workspace.status.objLoB == -Inf
                        workspace.status.absoluteGap = workspace.status.relativeGap = Inf
                    else
                        workspace.status.absoluteGap = workspace.status.objUpB - workspace.status.objLoB
                        workspace.status.relativeGap = workspace.status.absoluteGap/(1e-10 + abs(workspace.status.objUpB))
                    end
                end
            end


            # receive the number of solutions found and the global upper bound
            workspace.status.objUpB = workspace.sharedMemory.objectiveBounds[end]
            workspace.status.objLoB = min(workspace.status.objLoB,workspace.status.objUpB)
            workspace.status.numSolutions = workspace.sharedMemory.stats[1]

            # check if the new objective upper bound is lower that the user-defined cutoff
            if workspace.status.objUpB <= workspace.settings.objectiveCutoff
                # declare the cutoff inactive
                workspace.status.cutoffActive = false
            end

            # communicate the current local lower bound
            workspace.sharedMemory.objectiveBounds[processId] = workspace.status.objLoB

            # recompute optimality gaps
            if workspace.status.objUpB == Inf || workspace.status.objLoB == -Inf
                workspace.status.absoluteGap = workspace.status.relativeGap = Inf
            else
                workspace.status.absoluteGap = workspace.status.objUpB - workspace.status.objLoB
                workspace.status.relativeGap = workspace.status.absoluteGap/(1e-10 + abs(workspace.status.objUpB))
            end
        end
    end

    return
end
