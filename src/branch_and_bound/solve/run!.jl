# @Author: Massimo De Mauri <massimo>
# @Date:   2019-02-06T18:33:16+01:00
# @Email:  massimo.demauri@gmail.com
# @Project: OpenBB
# @Filename: run!.jl
# @Last modified by:   massimo
# @Last modified time: 2019-05-22T10:39:24+02:00
# @License: apache 2.0
# @Copyright: {{copyright}}


# insert a list of equivalent nodes into a queue
function insert_nodes!(queue::Array{BBnode,1},nodes::Array{BBnode,1},
                       priorityRule::Function,status::BBstatus;
                       unreliablePriority::Int=0)::Nothing

    if nodes[1].reliable || unreliablePriority == 0
    # normal queue insertion
        tmpIndex = 0
        for i in length(queue):-1:1
            if priorityRule(nodes[1],queue[i],status)
                tmpIndex = i
                break
            end
        end
        splice!(queue,tmpIndex+1:tmpIndex,nodes)

    elseif unreliablePriority == -1
        # put the new unreliable nodes at the bottom of the activeQueue to deal with them later
        splice!(workspace.activeQueue,1:0,out[2])

    elseif unreliablePriority == 1
        # put the new unreliable nodes at the top of the activeQueue to try to fastly get rid of them
        append!(workspace.activeQueue,out[2])
    else
        @error "wrong priority setting for unreliable nodes"
    end

    return
end





# branch and bound algorithm
function run!(workspace::BBworkspace)::Nothing

    # timing
    lastTimeCheckpoint = time()

    # update algorithm status
    workspace.status.description = "running"

    # set the algorithm in "active state"
    idle = false

    # main loop
    while !idle

        # update total time
        workspace.status.totalTime += time() - lastTimeCheckpoint; lastTimeCheckpoint = time()

        # stopping conditions
        if length(workspace.activeQueue) == 0 || # no more nodes in the queue
           workspace.status.totalTime >= workspace.settings.timeLimit || # time is up
           workspace.settings.custom_stopping_rule(workspace) || # custom stopping rule triggered
           workspace.status.absoluteGap <= workspace.settings.absoluteGapTolerance || # reached required absolute gap
           workspace.status.relativeGap <= workspace.settings.relativeGapTolerance || # reached required relative gap
           (workspace.settings.numSolutionsLimit > 0 && workspace.status.numSolutions >= workspace.settings.numSolutionsLimit) # the required number of solutions has been found

           # set the algorithm in idle state
           idle = true

       else # continue with branch and bound

           # update algorithm status and print it
           if workspace.settings.verbose
               if workspace.settings.iterationInfoFreq == 1 ||
                  mod(workspace.status.numRelaxationsSolved,workspace.settings.iterationInfoFreq) == 1
                   print_status(workspace)
              end
           end

           # pick a node to process from the activeQueue
           node = pop!(workspace.activeQueue)

           # pick a node to send to the neighbouring process from the activeQueue
           nodeToSend = NullBBnode()
           if workspace.globalInfo != nothing && # multiprocessing?
              !isready(workspace.outputChannel) && # send only one node per time
              length(workspace.activeQueue) > 0 && # there is a node to send
              workspace.activeQueue[end].objVal < workspace.status.objUpB
                nodeToSend = pop!(workspace.activeQueue)
                @async put!(workspace.outputChannel,nodeToSend)
           end

            if false # debug only
                println("node")
                println("- lBnd: ",node.branchLoBs)
                println("- uBnd: ",node.branchUpBs)
                println("- primal: ",node.primal)
                println("- objV: ",node.objVal," - avgF: ",node.avgFrac)
            end

            # solve the node
            out = solve_and_branch!(node,workspace)


            if false # debug only
                if out[1] == "infeasible"
                    println("\n  - infeasible")
                elseif out[1] == "solution"
                    println("\n  - solution")
                else
                    println("\n  - ",out[1])
                    for prb in out[2]
                        # println(" - prim: ",prb.primal[workspace.dscIndices])
                        println("   lBnd:",prb.branchLoBs)
                        println("   uBnd:",prb.branchUpBs)
                        println()
                    end
                end
            end


            if out[1] == "solution" && out[2][1].reliable # a reliable solution has been found
                # insert new solution into the solutionPool
                push!(workspace.solutionPool,out[2][1])

                # update the number of solutions found
                workspace.status.numSolutions += 1
                # update the objective upper bound
                workspace.status.objUpB = out[2][1].objVal

                # update the global objective upper bound and the number of solutions found
                if workspace.globalInfo != nothing && out[2][1].objVal < workspace.globalInfo[1]
                    workspace.globalInfo[1] = out[2][1].objVal
                    workspace.globalInfo[3] +=1
                end

            elseif out[1] == "solution"  # a not reliable solution has been found
                # store the obtained (not reliable) solution
                push!(workspace.unactivePool,out[2][1])

            elseif out[1] == "suboptimal" && workspace.settings.dynamicMode # in dynamic mode the suboptimal nodes cannot be completely eliminated
                # store the suboptimal node
                push!(workspace.unactivePool,out[2][1])

            elseif out[1] == "children"  # no solution found. two children nodes have been created

                # send one of the children to the neighbouring node
                if false && workspace.globalInfo != nothing && # multiprocessing?
                   !isready(workspace.outputChannel) && # send only one node per time
                   length(workspace.activeQueue) > 0 && # there is a node to send
                   workspace.activeQueue[end].objVal < workspace.status.objUpB

                     # send one of the children to the neighbouring node
                     put!(workspace.outputChannel,out[2][end])

                     # insert the other children in the activeQueue
                     insert_nodes!(workspace.activeQueue,out[2][1:end-1],workspace.settings.expansion_priority_rule,
                              workspace.status,unreliablePriority=workspace.settings.unreliable_subps_priority)

                else
                    # insert the new nodes into the activeQueue
                    insert_nodes!(workspace.activeQueue,out[2],workspace.settings.expansion_priority_rule,
                                  workspace.status,unreliablePriority=workspace.settings.unreliable_subps_priority)
                end

                # apply rounding heuristics
                if node.avgFrac <= workspace.settings.roundingHeuristicsThreshold
                    heuristicnode = simple_rounding_heuristics(node,workspace)
                    push!(workspace.activeQueue,heuristicnode)
                end
            end

            # recompute the lower bound if:
            # 1 - there is no lower bound
            # 2 - the node providing the lower bound has been removed from the activeQueue
            # 3 - the queue is empty
            if workspace.status.objLoB == -Inf ||
               node.objVal == workspace.status.objLoB ||
               (typeof(nodeToSend) != NullBBnode && nodeToSend.objVal == workspace.status.objLoB)

                newObjLoB = workspace.status.objUpB
                for i in length(workspace.activeQueue):-1:1
                    # if we have unreliable problems in the activeQueue we cannot update the lower bound
                    if !workspace.activeQueue[i].reliable
                        newObjLoB = workspace.status.objLoB
                        break
                    elseif workspace.activeQueue[i].objVal < newObjLoB
                        newObjLoB =  workspace.activeQueue[i].objVal
                    end
                end
                if workspace.status.objLoB > newObjLoB + workspace.settings.primalTolerance
                    println("branch and bound: the objective lower bound has decreased from "*string(workspace.status.objLoB)*" to "*string(newObjLoB)*"...")
                    println(workspace.status.objUpB)
                    println(length(workspace.activeQueue))
                end
                workspace.status.objLoB = newObjLoB
            end
        end

        # communicate with the other processes
        if workspace.globalInfo != nothing

            # check arrest conditions
            if idle && # nothing to do locally
               !isready(workspace.inputChannel) && # no nodes to pick
               workspace.globalInfo[2] == workspace.settings.numProcesses-1 # all the other workers are waiting. Start arrest procedure

                 # send a killer node to the neighbouring process
                 @async put!(workspace.outputChannel,KillerNode())

                 # go in waiting mode
                 workspace.globalInfo[2] += 1.

                 # wait for the arrival of a KillerNode
                 @assert take!(workspace.inputChannel) isa KillerNode

                 # exit waiting mode
                 workspace.globalInfo[2] -= 1.

            else

                # check if a new node is available
                while idle || isready(workspace.inputChannel)

                    if idle
                        # communicate the waiting state
                        workspace.globalInfo[2] += 1.
                    end

                    # take a new node from the input channel
                    newNode = take!(workspace.inputChannel)

                    if idle
                        # exit waiting state
                        workspace.globalInfo[2] -= 1.
                    end



                    if newNode isa KillerNode # oh no! a killer node! we have to stop!

                        # propagate the killerNode
                        @async put!(workspace.outputChannel,newNode)

                        # go idle
                        idle = true

                        # exit the inner loop
                        break

                    else # a normal node: insert it in the queue

                        # go active
                        idle = false

                        # insert the new node in the queue
                        insert_nodes!(workspace.activeQueue,[newNode],workspace.settings.expansion_priority_rule,
                                      workspace.status,unreliablePriority=workspace.settings.unreliable_subps_priority)

                        # update the objective lower bound
                        if newNode.objVal < workspace.status.objLoB
                            workspace.status.objLoB = newNode.objVal
                        end
                    end
                end
            end

            # update the number of solutions found and the upper bound
            workspace.status.objUpB = workspace.globalInfo[1]
            workspace.status.objLoB = min(workspace.status.objLoB,workspace.status.objUpB)
            workspace.status.numSolutions = workspace.globalInfo[3]
        end

        # recompute optimality gaps
        if workspace.status.objUpB == Inf || workspace.status.objLoB == -Inf
            workspace.status.absoluteGap = workspace.status.relativeGap = Inf
        else
            workspace.status.absoluteGap = workspace.status.objUpB - workspace.status.objLoB
            workspace.status.relativeGap = workspace.status.absoluteGap/abs(1e-10 + workspace.status.objUpB)
        end
    end

    ############################## termination ##############################
    if workspace.status.absoluteGap < workspace.settings.absoluteGapTolerance ||
       workspace.status.relativeGap < workspace.settings.relativeGapTolerance

        workspace.status.description = "optimalSolutionFound"
        if workspace.settings.verbose
            print_status(workspace)
            println(" Exit: Optimal Solution Found")
        end

    elseif length(workspace.activeQueue) == 0 && workspace.status.objUpB == Inf

        workspace.status.description = "infeasible"
        if workspace.settings.verbose
            print_status(workspace)
            println(" Exit: infeasibilty detected")
        end
    # elseif length(workspace.activeQueue) == 0 && length(workspace.solutionPool) > 0
    #
    #     workspace.status.description = "noReliableSolutionFound"
    #     if workspace.settings.verbose
    #         print_status(workspace)
    #         println(" Exit: no reliable solution found")
    #     end

    else
        workspace.status.description = "interrupted"
        if workspace.settings.verbose
            print_status(workspace)
            println(" Exit: interrupted")
        end
    end

    return
end
