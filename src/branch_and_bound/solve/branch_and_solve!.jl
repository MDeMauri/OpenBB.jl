# @Author: Massimo De Mauri <massimo>
# @Date:   2019-06-04T16:11:53+02:00
# @Email:  massimo.demauri@gmail.com
# @Filename: branch_and_solve!.jl
# @Last modified by:   massimo
# @Last modified time: 2019-11-21T18:40:39+01:00
# @License: apache 2.0
# @Copyright: {{copyright}}

function branch_and_solve!(node::BBnode,workspace::BBworkspace{T1,T2,T3})::Array{BBnode,1} where T1<:Problem where T2<:AbstractWorkspace where T3<:AbstractSharedMemory

    # store useful info on the node (before they get modified)
    nodeObjVal = node.objVal
    nodePrimal = copy(node.primal)

    # update the node to the latest version
    nodeJustUpdated = false
    if node.version < workspace.updatesRegister.masterVersion
        nodeJustUpdated = true
        if update!(node,workspace) == false # update! returns false if the node is discovered infeasible
            return BBnode[]
        end
    end

    # create a list of children
    if node.avgAbsFrac == 0.0 || isnan(node.objVal)
        children, branchIndices_dsc, presolveIndices = [node], [0], [0]
    else
        children, branchIndices_dsc, presolveIndices = branch!(node,workspace)
    end

    # solve all the children
    for k in 1:length(children)
        # Preprocess & solve
        if preprocess!(children[k],workspace,presolveIndices,
                       withBoundsPropagation=workspace.settings.withBoundsPropagation)
            solve_node!(children[k],workspace)
        else
            children[k].objVal = Inf
        end

        # update pseudoCosts
        if !nodeJustUpdated && branchIndices_dsc[k]>0 && children[k].reliable && children[k].objVal < Inf

            # compute objective and primal variation
            deltaObjective = max(children[k].objVal-nodeObjVal,workspace.settings.primalTolerance) # the max filters out small numerical errors
            deltaVariable = children[k].primal[workspace.problem.varSet.dscIndices[branchIndices_dsc[k]]] - nodePrimal[workspace.problem.varSet.dscIndices[branchIndices_dsc[k]]]

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


function branch!(node::BBnode,workspace::BBworkspace{T1,T2,T3})::Tuple{Array{BBnode,1},Array{Int,1}, Array{Int,1}} where T1<:Problem where T2<:AbstractWorkspace where T3<:AbstractSharedMemory

    # select a branching index
    branchIndex_dsc, priorityScores = branching_priority_rule(workspace.settings.branchingPriorityRule,
                                                       node.primal[workspace.problem.varSet.dscIndices],workspace.problem.varSet.pseudoCosts,
                                                       workspace.settings.primalTolerance)
    @assert branchIndex_dsc > 0

    # get the index of the branching variable wrt all the variables
    branchIndex = workspace.problem.varSet.dscIndices[branchIndex_dsc]

    # check if the selected variable belongs to a sos1 group
    sos1Branching = false
    if length(workspace.problem.varSet.sos1Groups)>0 && workspace.problem.varSet.sos1Groups[branchIndex_dsc] != 0

        # collect all the variables belonging to the same sos1 groups
        sos1Group = [i for i in 1:length(workspace.problem.varSet.dscIndices)
                          if workspace.problem.varSet.sos1Groups[i] == workspace.problem.varSet.sos1Groups[branchIndex_dsc] &&
                          node.varLoBs[workspace.problem.varSet.dscIndices[i]] != node.varUpBs[workspace.problem.varSet.dscIndices[i]]]

        sos1Branching = true
    end

    # actually branch
    if sos1Branching == true && length(sos1Group) > 1 # SOS1 branching

        # order the variables in the sos1 group by their priority score
        permute!(sos1Group,sortperm((@. abs(node.primal[workspace.problem.varSet.dscIndices[sos1Group]])),rev=true))

        # create a list of children nodes
        children = Array{BBnode}(undef,2)

        # first child
        children[1] = deepcopy(node)
        @. children[1].varLoBs[workspace.problem.varSet.dscIndices[sos1Group[1:2:end]]] = 0.
        @. children[1].varUpBs[workspace.problem.varSet.dscIndices[sos1Group[1:2:end]]] = 0.
        @. children[1].primal[workspace.problem.varSet.dscIndices[sos1Group[1:2:end]]] = 0.

        # second child
        children[2] = node
        @. children[2].varLoBs[workspace.problem.varSet.dscIndices[sos1Group[2:2:end]]] = 0.
        @. children[2].varUpBs[workspace.problem.varSet.dscIndices[sos1Group[2:2:end]]] = 0.
        @. children[2].primal[workspace.problem.varSet.dscIndices[sos1Group[2:2:end]]] = 0.

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


function solve_node!(node::BBnode,workspace::BBworkspace{T1,T2,T3})::Nothing where T1<:Problem where T2<:AbstractWorkspace where T3<:AbstractSharedMemory

    # solve the node
    # node status guide:
    # 0 -> solved
    # 1 -> infeasible
    # 2 -> unreliable
    # 3 -> error
    if any(@. node.varLoBs > node.varUpBs + workspace.settings.primalTolerance) ||
       any(@. node.cnsLoBs > node.cnsUpBs + workspace.settings.primalTolerance)
        ssStatus = 1
    else
        (ssStatus,info) = solve!(node,workspace.subsolverWS)
    end


    if ssStatus == 0
        node.pseudoObjective = node.objVal
        for (k,i) in enumerate(workspace.problem.varSet.dscIndices)
            node.pseudoObjective +=  min(workspace.problem.varSet.pseudoCosts[1][k,1]*(node.primal[i]-floor(node.primal[i]+workspace.settings.primalTolerance)),
                                         workspace.problem.varSet.pseudoCosts[1][k,2]*(ceil(node.primal[i]-workspace.settings.primalTolerance)-node.primal[i]))/length(workspace.problem.varSet.dscIndices)
        end
        node.reliable = true
    elseif ssStatus == 1
        node.pseudoObjective = Inf
        node.reliable = true
    elseif ssStatus == 2
        node.reliable = false
    else
        @error "Error in the subsolver"
    end

    # count how many subproblems we have solved
    workspace.status.numExploredNodes = workspace.status.numExploredNodes + 1

    # compute node average fractionality and pseudo_cost
    absoluteFractionality = @. abs(node.primal[workspace.problem.varSet.dscIndices] - round(node.primal[workspace.problem.varSet.dscIndices]))
    @. absoluteFractionality =  absoluteFractionality*(absoluteFractionality>workspace.settings.primalTolerance)
    node.avgAbsFrac = 2. * sum(absoluteFractionality)/length(workspace.problem.varSet.dscIndices)
    return
end
