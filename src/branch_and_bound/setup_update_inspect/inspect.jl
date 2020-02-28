
# @Author: Massimo De Mauri <massimo>
# @Date:   2019-02-25T17:07:48+01:00
# @Email:  massimo.demauri@gmail.com
# @Filename: inspect.jl
# @Last modified by:   massimo
# @Last modified time: 2019-09-25T17:21:24+02:00
# @License: LGPL-3.0
# @Copyright: {{copyright}}

# function to print BB status info to screen
function print_status(workspace::BBworkspace{T1,T2,T3})::Nothing where T1<:Problem where T2<:AbstractWorkspace where T3<:AbstractSharedMemory

    if workspace.sharedMemory isa NullSharedMemory

        globalObjectiveUpB = workspace.status.objUpB
        globalObjectiveLoB = workspace.status.objLoB
        globalAbsGap = workspace.status.absoluteGap
        globalRelGap = workspace.status.relativeGap

    else
        globalObjectiveUpB = workspace.sharedMemory.objectiveBounds[end]
        globalObjectiveLoB = minimum(workspace.sharedMemory.objectiveBounds[1:end-1])

        # recompute optimality gaps
        if globalObjectiveUpB == Inf || globalObjectiveLoB == -Inf
            globalAbsGap = globalRelGap =  Inf
        else
            globalAbsGap = globalObjectiveUpB - globalObjectiveLoB
            globalRelGap = globalAbsGap/(1e-10 + abs(globalObjectiveUpB))
        end

    end

    println("time: ",round(workspace.status.totalTime,digits = 2),
            " | best obj.: ", workspace.status.objUpB,
            " | best possible: ", globalObjectiveLoB,
            " | abs. gap.: ", globalAbsGap,
            " | rel. gap.: ",round(globalRelGap, digits = 2))

    return
end


# returns the best solution
function get_best_solution(workspace::BBworkspace{T1,T2,T3};localOnly::Bool=false)::AbstractBBnode where T1<:Problem where T2<:AbstractWorkspace where T3<:AbstractSharedMemory

    solution = NullBBnode()

    if length(workspace.solutionPool) != 0
        solution = workspace.solutionPool[end]
    end

    # check the other workers if needed/required
    if !(localOnly || workspace.sharedMemory isa NullSharedMemory) &&
       (solution isa NullBBnode || solution.objUpB > workspace.status.objLoB)
        for p in 2:workspace.settings.numProcesses
            node = remotecall_fetch(Main.eval,p,:(OpenBB.get_best_solution(workspace,localOnly=true)))
            if !(node isa NullBBnode) && node.objUpB == workspace.sharedMemory.objectiveBounds[end]
                solution = node
                break
            end
        end
    end
    return solution
end


# returns all the solutions
function get_all_solutions(workspace::BBworkspace{T1,T2,T3};localOnly::Bool=false)::Array{AbstractBBnode,1} where T1<:Problem where T2<:AbstractWorkspace where T3<:AbstractSharedMemory

    # collect all the local solutions
    solutions = copy(workspace.solutionPool)

    # check the other workers if needed/required
    if !(localOnly || workspace.sharedMemory isa NullSharedMemory)
        for p in 2:workspace.settings.numProcesses
            append!(solutions,remotecall_fetch(Main.eval,p,:(OpenBB.get_all_solutions(workspace,localOnly=true))))
        end
    end
    sort!(solutions,alg=MergeSort,rev=true,lt=(l,r)->lower_objective(l,r,workspace.status))
    return solutions
end


# returns the best node
function get_best_node(workspace::BBworkspace{T1,T2,T3};localOnly::Bool=false)::AbstractBBnode  where T1<:Problem where T2<:AbstractWorkspace where T3<:AbstractSharedMemory

    # define dummy best node
    bestNode = NullBBnode()

    if !localOnly && !(workspace.sharedMemory isa NullSharedMemory)

        nodes = Array{AbstractBBnode,1}(undef,workspace.settings.numProcesses)

        # call the local version of the function on the remote workers
        @sync for p in 2:workspace.settings.numProcesses
            @async nodes[p] = remotecall_fetch(Main.eval,p,:(OpenBB.get_best_node(workspace,localOnly=true)))
        end

        # call the local version of the function on the current process
        nodes[1] = get_best_node(workspace,localOnly=true)

        # choose the best of the returned nodes
        for node in nodes
            if bestNode isa NullBBnode || # there is no best node yet
               !(node isa NullBBnode) && ((node.avgFractionality==0 && bestNode.avgFractionality > 0) || # the new node is a solution while the best so far isn't
                                          expansion_priority_rule(workspace.settings.expansionPriorityRule,node,bestNode,workspace.status)) # the new node is better than the best so far

                # set the new node as the best one
                bestNode = node
           end
       end
   else
        # take the last solution generated if any
        if length(workspace.solutionPool) > 0
            bestNode = workspace.solutionPool[end]
        end
        # otherwise, take the first problem in the active queue
        if bestNode isa NullBBnode && length(workspace.activeQueue) > 0
            bestNode = workspace.activeQueue[end]
        end
        # otherwise, check the unactivePool for better nodes
        if bestNode isa NullBBnode && length(workspace.unactivePool) > 0
            for node in workspace.unactivePool
                if bestNode isa NullBBnode || expansion_priority_rule(workspace.settings.expansionPriorityRule,bestNode,node,workspace.status)
                    bestNode = node
                end
            end
        end
    end

    return deepcopy(bestNode)
end


# returns the number of variables
function get_numVariables(workspace::BBworkspace{T1,T2,T3})::Int where T1<:Problem where T2<:AbstractWorkspace where T3<:AbstractSharedMemory
    return get_size(workspace.problem.varSet)
end

# returns the number of discrete variables
function get_numDiscreteVariables(workspace::BBworkspace{T1,T2,T3})::Int where T1<:Problem where T2<:AbstractWorkspace where T3<:AbstractSharedMemory
    return get_numDiscrete(workspace.problem.varSet)
end

# returns the number of constraints
function get_numConstraints(workspace::BBworkspace{T1,T2,T3})::Int where T1<:Problem where T2<:AbstractWorkspace where T3<:AbstractSharedMemory
    return get_size(workspace.problem.cnsSet)
end

# ...
function get_constraints(workspace::BBworkspace{T1,T2,T3})::AbstractConstraintSet where T1<:Problem where T2<:AbstractWorkspace where T3<:AbstractSharedMemory
    return workspace.problem.cnsSet
end

# ...
function get_objective(workspace::BBworkspace{T1,T2,T3})::AbstractObjective where T1<:Problem where T2<:AbstractWorkspace where T3<:AbstractSharedMemory
    return workspace.problem.objFun
end

# this function returns the sparsity pattern of the constraint set
function get_constraints_sparsity(workspace::BBworkspace{T1,T2,T3})::Tuple{Array{Int,1},Array{Int,1}} where T1<:Problem where T2<:AbstractWorkspace where T3<:AbstractSharedMemory
    return get_sparsity(workspace.problem.cnsSet)
end

# this function returns the sparsity pattern a constraint in the constraint set
function get_constraint_sparsity(workspace::BBworkspace{T1,T2,T3},index::Int)::Array{Int,1}  where T1<:Problem where T2<:AbstractWorkspace where T3<:AbstractSharedMemory
    return get_sparsity(workspace.problem.cnsSet,index)
end

# this function returns the sparsity pattern of the objective function
function get_objective_sparsity(workspace::BBworkspace{T1,T2,T3})::Tuple  where T1<:Problem where T2<:AbstractWorkspace where T3<:AbstractSharedMemory
    return get_sparsity(workspace.problem.objFun)
end

#
function get_variableBounds(workspace::BBworkspace{T1,T2,T3})::Tuple{Array{Float64,1},Array{Float64,1}}  where T1<:Problem where T2<:AbstractWorkspace where T3<:AbstractSharedMemory
    return get_bounds(workspace.problem.varSet)
end

#
function get_constraintBounds(workspace::BBworkspace{T1,T2,T3})::Tuple{Array{Float64,1},Array{Float64,1}}  where T1<:Problem where T2<:AbstractWorkspace where T3<:AbstractSharedMemory
    return get_bounds(workspace.problem.cnsSet)
end

# return the status of the BB process
function get_status(workspace::BBworkspace{T1,T2,T3};localOnly::Bool=false)::BBstatus  where T1<:Problem where T2<:AbstractWorkspace where T3<:AbstractSharedMemory

    status = BBstatus(workspace.status)

    if !localOnly && !(workspace.sharedMemory isa NullSharedMemory)
        for p in 2:workspace.settings.numProcesses
            tmpStatus = remotecall_fetch(Main.eval,p,:(OpenBB.get_status(workspace,localOnly=true)))

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
        if status.objUpB == Inf || status.objLoB == -Inf
            status.absoluteGap = status.relativeGap = Inf
        else
            status.absoluteGap = status.objUpB - status.objLoB
            status.relativeGap = status.absoluteGap/(1e-10 + abs(status.objUpB))
        end
    end

    return status
end


# ...
function get_discreteIndices(workspace::T)::Array{Int,1} where T<:BBworkspace
    return workspace.problem.varSet.dscIndices
end

# ...
function get_sos1Groups(workspace::T)::Array{Int,1} where T<:BBworkspace
    return workspace.problem.varSet.sos1Groups
end

#...
function get_pseudoCosts(workspace::T)::Tuple{Array{Float64,2},Array{Int,2}} where T<:BBworkspace
    return workspace.problem.varSet.pseudoCosts
end
