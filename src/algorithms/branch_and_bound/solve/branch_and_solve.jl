# @Author: Massimo De Mauri <massimo>
# @Date:   2019-06-04T16:11:53+02:00
# @Email:  massimo.demauri@gmail.com
# @Filename: branch_and_solve!.jl
# @Last modified by:   massimo
# @Last modified time: 2021-02-11T21:33:00+01:00
# @License: apache 2.0
# @Copyright: {{copyright}}

function branch_and_solve!(node::BBnode,workspace::BBworkspace,nodeJustUpdated::Bool=false)::Vector{BBnode}


    # store useful info on the node (before they get modified)
    nodeObjInterval = (node.objLoB,node.objUpB)
    nodePrimal = copy(node.primal)

    # create a list of children
    if isempty(workspace.tree.active) && nodeObjInterval[1] == -Inf # root
        children, branchIndices_dsc, presolveIndices = [node], [0], [0]
    elseif nodeObjInterval[1] == -Inf # lower bound update failed
        children, branchIndices_dsc, presolveIndices = [node], [0], Int[]
    # elseif nodeObjInterval[2] == Inf # node to be re-solved
    #     children, branchIndices_dsc, presolveIndices = [node], [0], Int[]
    elseif node.avgFractionality == 0.0 # node already integer-feasible (to-resolve)
        children, branchIndices_dsc, presolveIndices = [node], [0], Int[]
    else # node to branch on
        children, branchIndices_dsc, presolveIndices = branch!(node,workspace)
    end


    # solve all the children
    for k in 1:length(children)
        # Preprocess & solve
        if isempty(presolveIndices) || preprocess_node!(children[k],workspace,presolveIndices,withBoundsPropagation=workspace.settings.withBoundsPropagation)
            solve_node!(children[k],workspace)
        else
            children[k].objLoB = children[k].objUpB = Inf
        end

        # update pseudoCosts
        if !nodeJustUpdated && nodeObjInterval[2] < Inf && branchIndices_dsc[k] > 0 && abs(children[k].objLoB) < Inf && children[k].reliable

            # compute objective and primal variation
            dscIndices = workspace.problem.varSet.dscIndices
            deltaObjective = max(children[k].objLoB-nodeObjInterval[1],workspace.settings.primalTolerance) # the max filters out small numerical errors
            deltaVariable = children[k].primal[dscIndices[branchIndices_dsc[k]]] - nodePrimal[dscIndices[branchIndices_dsc[k]]]

            if deltaVariable < -workspace.settings.primalTolerance

                # update the pseudoCost
                mu = 1/(workspace.problem.varSet.pseudoCosts[2][branchIndices_dsc[k],1]+1)
                workspace.problem.varSet.pseudoCosts[1][branchIndices_dsc[k],1] = (1-mu)*workspace.problem.varSet.pseudoCosts[1][branchIndices_dsc[k],1] - mu*deltaObjective/deltaVariable
                workspace.problem.varSet.pseudoCosts[2][branchIndices_dsc[k],1] += 1

            elseif deltaVariable > workspace.settings.primalTolerance

                # update the pseudoCost
                mu = 1/(workspace.problem.varSet.pseudoCosts[2][branchIndices_dsc[k],2]+1)
                workspace.problem.varSet.pseudoCosts[1][branchIndices_dsc[k],2] = (1-mu)*workspace.problem.varSet.pseudoCosts[1][branchIndices_dsc[k],2] + mu*deltaObjective/deltaVariable
                workspace.problem.varSet.pseudoCosts[2][branchIndices_dsc[k],2] += 1
            end
        end
    end

    # return the solved children
    return children
end


function branch!(node::BBnode,workspace::BBworkspace)::Tuple{Vector{BBnode},Vector{Int}, Vector{Int}}

    # collect the indices for the discrete variables
    dscIndices = workspace.problem.varSet.dscIndices

    # select a branching index
    branchIndex_dsc, bestScore = branching_priority_rule(workspace.settings.branchingPriorityRule,
                                                         node.primal[dscIndices],workspace.problem.varSet.pseudoCosts,
                                                         workspace.settings.primalTolerance)
    # the branching priority rule couldn't determine a
    # branching variable the node needs to be solved again
    # (this may happen when the bounds of the node have been set without solving it)
    if branchIndex_dsc == 0
        return [node], [0], Int[]
    end

    # get the index of the branching variable wrt all the variables
    branchIndex = dscIndices[branchIndex_dsc]

    # check if the selected variable belongs to a sos1 group
    sos1Branching = false
    if length(workspace.problem.varSet.sos1Groups)>0 && workspace.problem.varSet.sos1Groups[branchIndex_dsc] != 0

        # collect all the variables belonging to the same sos1 groups
        sos1Group = [i for i in 1:length(dscIndices)
                          if workspace.problem.varSet.sos1Groups[i] == workspace.problem.varSet.sos1Groups[branchIndex_dsc] &&
                          node.varLoBs[dscIndices[i]] != node.varUpBs[dscIndices[i]]]

        sos1Branching = true
    end

    # actually branch
    if sos1Branching == true && length(sos1Group) > 1 # SOS1 branching

        # order the variables in the sos1 group by their priority score
        sort!(sos1Group,lt=(l,r)->(abs(node.primal[dscIndices[l]])>=abs(node.primal[dscIndices[r]])))

        # create a list of children nodes
        children = Array{BBnode}(undef,2)

        # first child
        children[1] = deepcopy(node)
        @. children[1].varLoBs[dscIndices[sos1Group[1:2:end]]] = 0.
        @. children[1].varUpBs[dscIndices[sos1Group[1:2:end]]] = 0.
        @. children[1].primal[dscIndices[sos1Group[1:2:end]]] = 0.

        # second child
        children[2] = node
        @. children[2].varLoBs[dscIndices[sos1Group[2:2:end]]] = 0.
        @. children[2].varUpBs[dscIndices[sos1Group[2:2:end]]] = 0.
        @. children[2].primal[dscIndices[sos1Group[2:2:end]]] = 0.

        if length(sos1Group) == 3
            return children, [0,sos1Group[2]], sos1Group
        elseif length(sos1Group) == 2
            return children, [sos1Group[1],sos1Group[2]], sos1Group
        else
            return children, [0,0], sos1Group
        end

    else # standard branching

        # create a list of children nodes
        children = Array{BBnode}(undef,2)

        # first child
        children[1] = deepcopy(node)
        children[1].primal[branchIndex] = ceil(children[1].primal[branchIndex]-workspace.settings.primalTolerance)
        children[1].varLoBs[branchIndex] = children[1].primal[branchIndex]

        # second child
        children[2] = node
        children[2].primal[branchIndex] = floor(children[2].primal[branchIndex]+workspace.settings.primalTolerance)
        children[2].varUpBs[branchIndex] = children[2].primal[branchIndex]

        return children, [branchIndex_dsc,branchIndex_dsc], [branchIndex_dsc]

    end
end


function solve_node!(node::BBnode,workspace::BBworkspace)::Nothing

    # solve the node
    # node status guide:
    # 0 -> solved
    # 1 -> infeasible
    # 2 -> unreliable
    # 3 -> error
    if any(@. node.varLoBs > node.varUpBs + workspace.settings.primalTolerance) ||
       any(@. node.cnsLoBs > node.cnsUpBs + workspace.settings.primalTolerance)
        ssStatus = 1
        node.objLoB = node.objUpB = Inf
    else
        (ssStatus,info) = solve!(node,workspace.subsolverWS,objUpperLimit=get_suboptimality_threshold(workspace))
        workspace.status.numExploredNodes += 1
    end

    # round the discrete variables for ensuring correctness of branching
    dscIndices = workspace.problem.varSet.dscIndices
    @. node.primal[dscIndices] =
        clip(node.primal[dscIndices],
             node.varLoBs[dscIndices]-1e-1*workspace.settings.primalTolerance,
             node.varUpBs[dscIndices]+1e-1*workspace.settings.primalTolerance)

    if ssStatus == 0
        node.pseudoCost = 0.0
        for (k,i) in enumerate(dscIndices)
            node.pseudoCost += min(workspace.problem.varSet.pseudoCosts[1][k,1]*(node.primal[i]-floor(node.primal[i]+workspace.settings.primalTolerance)),
                                    workspace.problem.varSet.pseudoCosts[1][k,2]*(ceil(node.primal[i]-workspace.settings.primalTolerance)-node.primal[i]))
        end
        node.reliable = true
    elseif ssStatus == 1
        node.reliable = true
    elseif ssStatus == 2
        node.reliable = false
    else
        @error "Error in the subsolver"
    end

    # compute node average fractionality
    absoluteFractionality = @. abs(node.primal[dscIndices] - round(node.primal[dscIndices]))
    @. absoluteFractionality =  absoluteFractionality*(absoluteFractionality>workspace.settings.primalTolerance)
    node.avgFractionality = 2.0*sum(absoluteFractionality)/Float(length(absoluteFractionality))
    return
end
