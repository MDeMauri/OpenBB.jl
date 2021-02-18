# @Author: Massimo De Mauri <massimo>
# @Date:   2019-02-06T18:33:16+01:00
# @Email:  massimo.demauri@gmail.com
# @Project: OpenBB
# @Filename: run!.jl
# @Last modified by:   massimo
# @Last modified time: 2021-02-12T19:25:44+01:00
# @License: LGPL-3.0
# @Copyright: {{copyright}}


# branch and bound algorithm
function run!(workspace::BBworkspace)::Nothing
    # timing
    lastTimeCheckpoint = time()

    # id of the current process
    processId = myid()

    # apply preprocessing on the root
    if workspace.status.description == "new" && processId == 1
        root = workspace.tree.active[end][1]
        if !update!(root,workspace)[2] || !preprocess_root!(root,workspace)
            @warn "Root problem found infeasible"
            workspace.status.objUpB = workspace.status.objLoB = Inf
            # update elapsed time
            workspace.status.totalTime += time() - lastTimeCheckpoint
            return
        end
    end

    # update algorithm status
    workspace.status.description = "running"

    # set the algorithm in "active state"
    idle = false

    # print the initial algorithm status
    printCountdown = 0.0

    # allowed time for inter-process communication
    communicationTimeout = 1e-2

    # keep in memory how many nodes the algorithm has explored since the last time
    # it has sent a node to the neighbouring process
    nodeSharingCountdown = workspace.settings.nodeSharingPeriod-1

    # keep in memory the objective of the node in the outputChannel
    objNodeInOutputChannel = Inf

    # communicate the initial objective lower-bound to the other nodes
    if !(workspace.sharedMemory isa NullSharedMemory)
        # communicate the current local lower bound
        workspace.sharedMemory.localObjLoBs[processId] = workspace.status.objLoB
    end

    # main loop
    while true

        # get elapsed time for the last iteration
        iterationTime = time() - lastTimeCheckpoint
        lastTimeCheckpoint = time()
        # update total time
        workspace.status.totalTime += iterationTime
        printCountdown -= iterationTime

        # communicate wtith other processes
        if !(workspace.sharedMemory isa NullSharedMemory)
            if isready(workspace.sharedMemory.inputChannel)

                # take a new node from the input channel
                newNode = take!(workspace.sharedMemory.inputChannel;timeout=communicationTimeout)

                # if the new node is a branch and bound node insert it in the tree
                if newNode isa BBnode
                    if idle && !is_suboptimal(newNode,workspace)
                        # insert node on top of the active queue
                        push_node!(workspace,:active,newNode)
                        # do not send this node to the neighbour
                        nodeSharingCountdown = max(nodeSharingCountdown,1)
                        # go active
                        idle = false
                        # recompute lower bound
                        workspace.status.objLoB = min(workspace.status.objLoB,newNode.objLoB)
                    else
                        # insert normally
                        insert_node!(workspace,newNode)
                    end

                elseif newNode isa KillerNode
                    if newNode.code == 0
                        # declare the error in the status
                        workspace.status.somethingWrong = true
                        workspace.status.description = "error"
                    end
                    # propagate the killer node
                    put!(workspace.sharedMemory.outputChannel,newNode,mode=:forced)
                    # terminate the process going arrestable
                    workspace.sharedMemory.arrestable[processId] = true
                    # update elapsed time
                    workspace.status.totalTime += time() - lastTimeCheckpoint
                    return
                end

                # check if there is still a node in the outputChannel
                if isempty(workspace.sharedMemory.outputChannel)
                    objNodeInOutputChannel = Inf
                end

            end

            # receive the number of solutions found and the global upper bound
            workspace.status.objUpB = workspace.sharedMemory.globalObjUpB[1][1]
            workspace.status.numSolutions = workspace.sharedMemory.stats[1]

            # communicate the current local lower bound
            workspace.sharedMemory.localObjLoBs[processId] = min(workspace.status.objLoB,objNodeInOutputChannel)

        end

        # recompute optimality gaps
        (workspace.status.absoluteGap,workspace.status.relativeGap) =
            compute_gaps(min(workspace.status.objLoB,objNodeInOutputChannel),workspace.status.objUpB,workspace.settings.gapsComputationMode)

        # print algorithm status
        if printCountdown <= 0
            if workspace.settings.verbose && processId == 1
                print_status(workspace)
            elseif processId == 1
                # I do not know why it is necessary, but it is:
                # otherwise the other processes do not work properly
                print("")
            end
            printCountdown = workspace.settings.statusInfoPeriod
        end

        # ideling conditions
        if !idle # we are not already idle
            if length(workspace.tree.active) == 0 || # no more nodes in the queue
               workspace.status.objLoB >= get_suboptimality_threshold(workspace) || # there exist no other solution within the cutoff or we satisfy the gap tolerances
               workspace.status.totalTime >= workspace.settings.timeLimit || # time is up
               workspace.settings.customInterruptionRule(workspace) || # custom stopping rule triggered
               (workspace.settings.numSolutionsLimit > 0 && workspace.status.numSolutions >= workspace.settings.numSolutionsLimit) # the required number of solutions has been found

                # set the algorithm in idle state
                idle = true
            end
        end

        # perform branch and bound iteration
        if !idle

            # pick a node to process from the tree.active
            (node,~) = pop!(workspace.tree.active)

            # store useful info for node and workspace
            nodeObjLoB = node.objLoB
            oldObjUpB = workspace.status.objUpB
            oldNumSolutions = workspace.status.numSolutions

            # update the node if necessary
            (nodeJustUpdated,nodeFeasible) = update!(node,workspace)
            if !nodeFeasible
                # round the discrete variables for ensuring correctness of future branching
                dscIndices_ = workspace.problem.varSet.dscIndices
                @. node.primal[dscIndices_] = clip(node.primal[dscIndices_],node.varLoBs[dscIndices_],node.varUpBs[dscIndices_])
                # mark the node as infeasible
                node.objLoB = node.objUpB = Inf
            end

            if is_suboptimal(node,workspace)

                if node.objLoB > workspace.settings.objectiveCutoff # the node violates the cutoff
                    # declare the cutoff active
                    workspace.status.cutoffActive = true
                end

                if workspace.settings.conservativismLevel >= 1 && !node.heuristic
                    # store the node as inactive
                    insert_node!(workspace,:suboptimals,node)
                end

            else # perform B&B iteration

                ######## Heuristic step ########
                # check existence of heuristics, feasibility/fractionality of the node and triggering condition
                if !(workspace.heuristicsWS isa NullBBheuristicsWorkspace) &&
                    nodeFeasible && !nodeJustUpdated && node.avgFractionality > 0.0 &&
                    check_heuristicsTriggerCondition(node,workspace)

                    # apply heuristics
                    hNodes = get_heuristic_nodes(node,workspace.heuristicsWS,workspace)

                    for hNode in hNodes
                        solve_node!(hNode,workspace)
                        insert_node!(workspace,hNode)
                    end
                end


                ####### Normal Step ########
                # reduce the nodeSharingCountdown
                nodeSharingCountdown -= 1

                # create a list of children nodes
                if nodeFeasible
                    children = branch_and_solve!(node,workspace,nodeJustUpdated)

                    # if we have more than one child and we are multiprocessing,
                    # send the second child to the next process (if possible).
                    if !(workspace.sharedMemory isa NullSharedMemory) && # we are multiprocessing
                       length(children) > 1                           && # we have more than one child
                       isempty(workspace.sharedMemory.outputChannel)  && # the output channel is free
                       !islocked(workspace.sharedMemory.outputChannel)   # the output channel is unlocked

                        # send the new node to the neighbouring process
                        put!(workspace.sharedMemory.outputChannel,children[2];timeout=communicationTimeout)

                        # remember the objective of the node in the output channel
                        objNodeInOutputChannel = node.objLoB

                        # remove the child from the list of children
                        deleteat!(children,2)
                    end

                    # insert the children in the BBtree (checking for solutions and suboptimals)
                    for child in children insert_node!(workspace,child) end

                else
                    insert_node!(workspace,node)
                end
            end

            # communicate the new objective upper-bound if necessary
            # update the global objective upper bound and the number of solutions found
            if !(workspace.sharedMemory isa NullSharedMemory)
                if oldNumSolutions != workspace.status.numSolutions
                    workspace.sharedMemory.stats[1] += workspace.status.numSolutions - oldNumSolutions
                end
                if oldObjUpB != workspace.status.objUpB && workspace.status.objUpB < workspace.sharedMemory.globalObjUpB[1][1]
                    workspace.sharedMemory.globalObjUpB[1] = (workspace.status.objUpB,Int8(myid()))
                end
            end


            # recompute the objective lower bound if needed
            if nodeObjLoB == workspace.status.objLoB || workspace.status.objLoB == -Inf
                if workspace.status.cutoffActive
                    newObjLoB = min(workspace.settings.objectiveCutoff,workspace.status.objUpB)
                else
                    newObjLoB = workspace.status.objUpB
                end
                for i in length(workspace.tree.active):-1:1
                    # if we have unreliable problems in the tree.active we cannot update the lower bound
                    if !workspace.tree.active[i][1].reliable
                        newObjLoB = workspace.status.objLoB
                        # update elapsed time
                        workspace.status.totalTime += time() - lastTimeCheckpoint
                        break
                    elseif newObjLoB > workspace.tree.active[i][1].objLoB
                        newObjLoB = workspace.tree.active[i][1].objLoB
                    end
                end

                # warn the user if the objective lower bound has significantly decreased (something might be wrong)
                if workspace.status.objLoB > newObjLoB + workspace.settings.primalTolerance
                    @warn "Branch and Bound: the objective lower bound has decreased of "*string(workspace.status.objLoB-newObjLoB)
                end

                # update the objective lower-bound
                workspace.status.objLoB = newObjLoB
            end

        # in single processing we stop as soon as the algorithm goes idle
        elseif workspace.sharedMemory isa NullSharedMemory
            # update elapsed time
            workspace.status.totalTime += time() - lastTimeCheckpoint
            return

        # in multi processing, if idle, the algorithm starts to wait for
        # for others to go idle or for a node to be sent from the neighbour,
        # and also for the last node it sent to be picked up
        elseif isempty(workspace.sharedMemory.inputChannel)

            # keep track of the time spent in waiting
            waitingStartTime = time()
            workspace.status.description = "waiting"
            nodeInOutput = true
            hasToStop = false
            while isempty(workspace.sharedMemory.inputChannel)
                if nodeInOutput && isempty(workspace.sharedMemory.outputChannel)
                    nodeInOutput = false
                    workspace.sharedMemory.arrestable[processId] = true
                elseif all(workspace.sharedMemory.arrestable)
                    workspace.status.description = "done"
                    hasToStop = true
                    # update elapsed time
                    workspace.status.totalTime += time() - lastTimeCheckpoint
                    return
                end
                sleep(0.001)
            end
            # undeclare the node ready to stop
            if !nodeInOutput
                workspace.sharedMemory.arrestable[processId] = false
            end
            workspace.status.description = "running"
            # keep track of the time spent waiting and of the total time
            workspace.status.waitingTime += time() - waitingStartTime
        end
    end
end
