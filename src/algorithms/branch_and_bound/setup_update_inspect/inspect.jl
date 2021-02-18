
# @Author: Massimo De Mauri <massimo>
# @Date:   2019-02-25T17:07:48+01:00
# @Email:  massimo.demauri@gmail.com
# @Filename: inspect.jl
# @Last modified by:   massimo
# @Last modified time: 2021-02-12T14:45:16+01:00
# @License: LGPL-3.0
# @Copyright: {{copyright}}

# function to print BB status info to screen
function print_status(workspace::BBworkspace)::Nothing
    if workspace.sharedMemory isa NullSharedMemory
        globalObjUpB = workspace.status.objUpB
        globalObjLoB = workspace.status.objLoB
        globalAbsGap = workspace.status.absoluteGap
        globalRelGap = workspace.status.relativeGap
    else
        globalObjUpB = workspace.sharedMemory.globalObjUpB[1][1]
        globalObjLoB = minimum(workspace.sharedMemory.localObjLoBs)
        (globalAbsGap,globalRelGap) = compute_gaps(globalObjLoB,globalObjUpB,workspace.settings.gapsComputationMode)
    end

    println(" - time: ",round(workspace.status.totalTime,digits = 2),
            " | best obj.: ", globalObjUpB,
            " | best possible: ", globalObjLoB,
            " | abs. gap.: ", globalAbsGap,
            " | rel. gap.: ",round(globalRelGap, digits = 2))

    return
end

function get_numFeasibleNodes(workspace::BBworkspace;localOnly::Bool=false)::Int
	if !localOnly && !(workspace.sharedMemory isa NullSharedMemory)
		counts = Vector{Int}(undef,workspace.settings.numProcesses)
		counts[1] = get_numFeasibleNodes(workspace,localOnly=true)
		@sync for (k,id) in enumerate(workers()[1:workspace.settings.numProcesses-1])
			counts[k+1] = remotecall_fetch(OpenBB.eval,id,:(get_numFeasibleNodes(workspace,localOnly=true)))
		end
		return sum(counts)
	else
		return length(workspace.tree.solutions)
	end
end

function get_numActiveNodes(workspace::BBworkspace;localOnly::Bool=false)::Int
	if !localOnly && !(workspace.sharedMemory isa NullSharedMemory)
		counts = Vector{Int}(undef,workspace.settings.numProcesses)
		counts[1] = get_numActiveNodes(workspace,localOnly=true)
		@sync for (k,id) in enumerate(workers()[1:workspace.settings.numProcesses-1])
			counts[k+1] = remotecall_fetch(OpenBB.eval,id,:(get_numActiveNodes(workspace,localOnly=true)))
		end
		return sum(counts)
	else
		return length(workspace.tree.active)
	end
end


function get_numInactiveNodes(workspace::BBworkspace;localOnly::Bool=false)::Int
	if !localOnly && !(workspace.sharedMemory isa NullSharedMemory)
		counts = Vector{Int}(undef,workspace.settings.numProcesses)
		counts[1] = get_numInactiveNodes(workspace,localOnly=true)
		@sync for (k,id) in enumerate(workers()[1:workspace.settings.numProcesses-1])
			counts[k+1] = remotecall_fetch(OpenBB.eval,id,:(get_numInactiveNodes(workspace,localOnly=true)))
		end
		return sum(counts)
	else
		return length(workspace.tree.suboptimals) + length(workspace.tree.infeasibles)
	end
end


# returns the best solution node
function get_best_feasible_node(workspace::BBworkspace;localOnly::Bool=false)::AbstractBBnode
    if !localOnly && !(workspace.sharedMemory isa NullSharedMemory)
        if workspace.sharedMemory.globalObjUpB[1][2] == 0
            return NullBBnode()
        elseif workspace.sharedMemory.globalObjUpB[1][2] == myid()
            return get_best_feasible_node(workspace,localOnly=true)
        else
            return remotecall_fetch(OpenBB.eval,workspace.sharedMemory.globalObjUpB[1][2],:(OpenBB.get_best_feasible_node(workspace,localOnly=true)))
        end
    else
        if !isempty(workspace.tree.solutions)
            return deepcopy(workspace.tree.solutions[end][1])
        else
            return NullBBnode()
        end
    end
end


# returns all the solutions
function get_all_feasible_nodes(workspace::BBworkspace;localOnly::Bool=false)::Vector{BBnode}
    # collect all the local solutions
    solutions = @. deepcopy(getindex(workspace.tree.solutions,1))

    # check the other workers if needed/required
    if !(localOnly || workspace.sharedMemory isa NullSharedMemory)
        for p in workers()[1:workspace.settings.numProcesses-1]
            append!(solutions,remotecall_fetch(OpenBB.eval,p,:(OpenBB.get_all_feasible_nodes(workspace,localOnly=true))))
        end
    end
    collect_score_and_sort!(solutions,node->node.objLoB,algorithm=MergeSort,reverse=true)
    return solutions
end


# returns the best node
function get_best_node(workspace::BBworkspace;localOnly::Bool=false)::BBnode

    # if we have a solution return it
    if workspace.status.objUpB < Inf
        return get_best_feasible_node(workspace,localOnly=localOnly)
    end

    # otherwise, look for general nodes
    if !localOnly && !(workspace.sharedMemory isa NullSharedMemory)

        candidates = Vector{AbstractBBnode}(undef,workspace.settings.numProcesses)

        # call the local version of the function on the remote workers
        @sync for (k,p) in enumerate(workers()[1:workspace.settings.numProcesses-1])
            @async candidates[k] = remotecall_fetch(OpenBB.eval,p,:(OpenBB.get_best_node(workspace,localOnly=true)))
        end

        # call the local version of the function on the current process
        candidates[end] = get_best_node(workspace,localOnly=true)

        # choose the best of the returned nodes
        filter!(x->!(x isa NullBBnode),candidates)
        if !isempty(candidates)
			bestNode = candidates[1]
			bestPriority = expansion_priority(candidates[1],workspace)
			for node in candidates[2:end]
				priority_ = expansion_priority(node,workspace)
				if priority_ > bestPriority
					bestNode = node
					bestPriority = priority_
				end
			end
		else
            bestNode = NullBBnode()
        end
   else
        if !isempty(workspace.tree.active) # Take the first problem in the active queue
            bestNode = workspace.tree.active[end][1]
        else # otherwise, return a null node
            bestNode = NullBBnode()
        end
    end

    return deepcopy(bestNode)
end


# returns the number of variables
function get_numVariables(workspace::BBworkspace)::Int
    return get_size(workspace.problem.varSet)
end

# returns the number of discrete variables
function get_numDiscrete(workspace::BBworkspace)::Int
    return get_numDiscrete(workspace.problem.varSet)
end

# returns the number of constraints
function get_numConstraints(workspace::BBworkspace)::Int
    return get_size(workspace.problem.cnsSet)
end

# ...
function get_constraints(workspace::BBworkspace)::ConstraintSet
    return workspace.problem.cnsSet
end

# ...
function get_objective(workspace::BBworkspace)::ObjectiveFunction
    return workspace.problem.objFun
end

# this function returns the list of variables occurring in the constraint set
function get_constraints_dependency(workspace::BBworkspace)::Vector{Vector{Int}}
    return get_dependency(workspace.problem.cnsSet)
end


# this function returns the list of variables occurring in a constraint in the constraint set
function get_constraint_dependency(workspace::BBworkspace,index::Int)::Vector{Int}
    return get_dependency(workspace.problem.cnsSet,index)
end

# this function returns the list of variables occurring in objective function
function get_objective_dependency(workspace::BBworkspace)::Vector{Int}
    return get_dependency(workspace.problem.objFun)
end

#
function get_variableBounds(workspace::BBworkspace)::Tuple{Vector{Float},Vector{Float}}
    return get_bounds(workspace.problem.varSet)
end

#
function get_constraintBounds(workspace::BBworkspace)::Tuple{Vector{Float},Vector{Float}}
    return get_bounds(workspace.problem.cnsSet)
end

# return the status of the BB process
function get_status(workspace::BBworkspace;localOnly::Bool=false)::BBstatus
    status = deepcopy(workspace.status)

    if !localOnly && !(workspace.sharedMemory isa NullSharedMemory)
        for p in workers()[1:workspace.settings.numProcesses-1]
            tmpStatus = remotecall_fetch(OpenBB.eval,p,:(OpenBB.get_status(workspace,localOnly=true)))

            # collect the objective bounds
            status.objLoB = min(status.objLoB,tmpStatus.objLoB)
            status.objUpB = min(status.objUpB,tmpStatus.objUpB)

            # collect reliability status
            status.reliable = status.reliable && tmpStatus.reliable

            # collect timings
            status.totalTime = max(status.totalTime,tmpStatus.totalTime)
            status.waitingTime = max(status.waitingTime,tmpStatus.waitingTime)

            # collect results statistics
            status.numSolutions += tmpStatus.numSolutions
            status.numExploredNodes += tmpStatus.numExploredNodes
        end

        # re-compute optimality gaps
        (status.absoluteGap,status.relativeGap) = compute_gaps(status.objLoB,status.objUpB,workspace.settings.gapsComputationMode)
    end

    return status
end


# ...
function get_discreteIndices(workspace::BBworkspace)::Vector{Int}
    return workspace.problem.varSet.dscIndices
end

# ...
function get_sos1Groups(workspace::BBworkspace)::Vector{Int}
    return workspace.problem.varSet.sos1Groups
end

#...
function get_pseudoCosts(workspace::BBworkspace)::Tuple{Matrix{Float},Matrix{Int}}
    return workspace.problem.varSet.pseudoCosts
end


# returns a decision on the suboptimality of a node depending on the status of the algorithm
function is_suboptimal(node::BBnode,workspace::BBworkspace)::Bool
    return node.objLoB >= get_suboptimality_threshold(workspace::BBworkspace)
end

# returns objective value after which a node is considered suboptimal
function get_suboptimality_threshold(workspace::BBworkspace)::Float
	if workspace.status.objUpB < Inf
	    threshold = compute_suboptimality_threshold(workspace.status.objUpB,
	                                                workspace.settings.absoluteGapTolerance,
	                                                workspace.settings.relativeGapTolerance,
	                                                workspace.settings.gapsComputationMode)
		return min(threshold,workspace.settings.objectiveCutoff)
	else
		return workspace.settings.objectiveCutoff
	end
end
