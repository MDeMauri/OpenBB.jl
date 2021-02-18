# @Author: Massimo De Mauri <massimo>
# @Date:   2019-08-13T17:42:17+02:00
# @Email:  massimo.demauri@gmail.com
# @Filename: insert_node.jl
# @Last modified by:   massimo
# @Last modified time: 2021-02-09T20:37:39+01:00
# @License: LGPL-3.0
# @Copyright: {{copyright}}



# insert a node into a node container
function insert_node!(workspace::BBworkspace,poolSymbol::Symbol,node::BBnode)::Nothing

    nodePool = getfield(workspace.tree,poolSymbol)

    if node.reliable || workspace.settings.unreliablesPriority == 0
    # normal queue insertion
        insertionPoint = 1
        nodeTuple = (node,expansion_priority(node,workspace))
        for i in length(nodePool):-1:1
            if nodeTuple[2] >= nodePool[i][2]
                insertionPoint = i+1
                break
            end
        end
        insert!(nodePool,insertionPoint,nodeTuple)

    elseif workspace.settings.unreliablesPriority == -1
        # put the new unreliable nodes at the bottom of the activeQueue to deal with them later
        insert!(nodePool,0,(node,-Inf))

    elseif workspace.settings.unreliablesPriority == 1
        # put the new unreliable nodes at the top of the activeQueue to try to fastly get rid of them
        append!(nodePool,(node,Inf))
    else
        @error "wrong priority setting for unreliable nodes"
    end

    return
end


function push_node!(workspace::BBworkspace,poolSymbol::Symbol,node::BBnode)::Nothing
    nodeTuple = (node,expansion_priority(node,workspace))
    nodePool = push!(getfield(workspace.tree,poolSymbol),nodeTuple)
    return
end

# insert a node into in one of the node containers
function insert_node!(workspace::BBworkspace,node::BBnode)::Nothing

    if node.objLoB == -Inf # not enough info on the node, reinsert it into active
        insert_node!(workspace,:active,node)
        workspace.status.objLoB = -Inf

    elseif node.objLoB == Inf # infeasible node

        if workspace.settings.conservativismLevel == 2 && !node.heuristic
            # store the node as inactive
            insert_node!(workspace,:infeasibles,node)
        end

    elseif node.avgFractionality == 0.0 &&  # a new solution has been found
           node.objUpB < workspace.settings.objectiveCutoff &&
           node.objUpB < workspace.status.objUpB

        if is_mixedBinary(workspace.problem.varSet) && lookup(workspace.blackList,node.primal[workspace.problem.varSet.dscIndices]) # the solution is blacklisted

            # collect info
            dscIndices = get_discreteIndices(workspace.problem)
            blackListedAssignment = round.(node.primal[dscIndices])


            # remove the blacklisted solution from the feasible set of the current node
            A = spzeros(1,get_size(workspace.problem.varSet))
            @. A[1,dscIndices] = 2.0*blackListedAssignment - 1.0
            add_cuts!(node,A=A,loBs=[-Inf],upBs=[sum(blackListedAssignment)-1.0],forBlackList=true)

            # reinsert the node into active
            node.objUpB = Inf
            node.avgFractionality = 1.0
            insert_node!(workspace,:active,node)

            # generate a new node representing only the blacklisted solution
            blNode = deepcopy(node)
            dscIndices = get_discreteIndices(workspace.problem.varSet)
            @. blNode.varLoBs[dscIndices] = blNode.varUpBs[dscIndices] = round(blNode.primal[dscIndices])

            # put the node in the blacklisted pool
            push!(workspace.tree.blacklisted,blNode)

            # declare the blackList active
            workspace.status.blackListActive = true

        elseif node.reliable || workspace.settings.acceptUnreliableSolutions # the solution is reliable or the algorithm is set to accept unreliable solutions

            # insert new solution into the solutionPool
            push_node!(workspace,:solutions,node)

            # update the number of solutions found
            workspace.status.numSolutions += 1

            # update the objective upper bound
            workspace.status.objUpB = node.objUpB

        else # the solution is not reliable

            # the unreliability of the node cannot be solved... duly noted!
            workspace.status.reliable = false

            # store the unreliable solution in the inactivePool
            insert_node!(workspace,:suboptimals,node)
        end

    elseif is_suboptimal(node,workspace) # suboptimal node

        if node.objLoB > workspace.settings.objectiveCutoff
            # declare the cutoff active
            workspace.status.cutoffActive = true
        end

        if workspace.settings.conservativismLevel >= 1 && !node.heuristic
            # store the node as inactive
            insert_node!(workspace,:suboptimals,node)
        end

    else # default behaviour
        insert_node!(workspace,:active,node)
        # update the objective lower bound
        workspace.status.objLoB = min(workspace.status.objLoB,node.objLoB)
    end

    return
end
