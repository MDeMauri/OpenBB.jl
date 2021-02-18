# @Author: Massimo De Mauri <massimo>
# @Date:   2020-11-26T11:34:22+01:00
# @Email:  massimo.demauri@gmail.com
# @Filename: inspect.jl
# @Last modified by:   massimo
# @Last modified time: 2021-01-13T16:10:49+01:00
# @License: LGPL-3.0
# @Copyright: {{copyright}}


# function to print HBB status info to screen
function print_status(workspace::HBBworkspace)::Nothing

    globalObjUpB = workspace.sharedMemory.globalObjUpB[1][1]
    println("time: ",round(workspace.status.totalTime,digits = 2),
            " | best obj.: ",workspace.status.objUpB,
            " | best possible: ",workspace.status.objLoB,
            " | abs. gap.: ",workspace.status.absoluteGap,
            " | rel. gap.: ",round(workspace.status.relativeGap,digits = 2))

    return
end

function get_numFeasibleNodes(workspace::HBBworkspace)::Int
    return length(workspace.feasibleNodes)
end

function get_best_feasible_node(workspace::HBBworkspace)::AbstractBBnode
    if isempty(workspace.feasibleNodes)
        return NullBBnode()
    else
        return workspace.feasibleNodes[1]
    end
end

function get_all_feasible_nodes(workspace::HBBworkspace)::Vector{BBnode}
    return deepcopy(workspace.feasibleNodes)
end


# returns the number of variables
function get_numVariables(workspace::HBBworkspace)::Int
    return get_size(workspace.problem.varSet)
end

# returns the number of discrete variables
function get_numDiscrete(workspace::HBBworkspace)::Int
    return get_numDiscrete(workspace.problem.varSet)
end

# returns the number of constraints
function get_numConstraints(workspace::HBBworkspace)::Int
    return get_size(workspace.problem.cnsSet)
end

# ...
function get_constraints(workspace::HBBworkspace)::ConstraintSet
    return workspace.problem.cnsSet
end

# ...
function get_objective(workspace::HBBworkspace)::ObjectiveFunction
    return workspace.problem.objFun
end

# this function returns the list of variables occurring in the constraint set
function get_constraints_dependency(workspace::HBBworkspace)::Vector{Vector{Int}}
    return get_dependency(workspace.problem.cnsSet)
end


# this function returns the list of variables occurring in a constraint in the constraint set
function get_constraint_dependency(workspace::HBBworkspace,index::Int)::Vector{Int}
    return get_dependency(workspace.problem.cnsSet,index)
end

# this function returns the list of variables occurring in objective function
function get_objective_dependency(workspace::HBBworkspace)::Vector{Int}
    return get_dependency(workspace.problem.objFun)
end

#
function get_variableBounds(workspace::HBBworkspace)::Tuple{Vector{Float},Vector{Float}}
    return get_bounds(workspace.problem.varSet)
end

#
function get_constraintBounds(workspace::HBBworkspace)::Tuple{Vector{Float},Vector{Float}}
    return get_bounds(workspace.problem.cnsSet)
end


# ...
function get_discreteIndices(workspace::HBBworkspace)::Vector{Int}
    return workspace.problem.varSet.dscIndices
end

# ...
function get_sos1Groups(workspace::HBBworkspace)::Vector{Int}
    return workspace.problem.varSet.sos1Groups
end

function get_status(workspace::HBBworkspace;localOnly::Bool=false)::HBBstatus
    return deepcopy(workspace.status)
end


# returns a decision on the suboptimality of a node depending on the status of the algorithm
function is_suboptimal(node::BBnode,workspace::HBBworkspace)::Bool

    if node.objLoB >= workspace.settings.objectiveCutoff - workspace.settings.primalTolerance
        return true
    elseif workspace.status.objUpB < Inf && node.objLoB > -Inf
        # compute the gaps assuming that the node is the only left in the tree,
        # and check if we would stop BB in that condition.
        (absGap_,relGap_) = compute_gaps(node.objLoB,workspace.status.objUpB,workspace.settings.gapsComputationMode)
        if absGap_ <= workspace.settings.absoluteGapTolerance || relGap_ <= workspace.settings.relativeGapTolerance
            return true
        else
            return false
        end
    else
        return false
    end
end

# returns objective value after which a node is considered suboptimal
function get_suboptimality_threshold(workspace::HBBworkspace)::Float
    threshold = compute_suboptimality_threshold(workspace.status.objUpB,
                                                workspace.settings.absoluteGapTolerance,
                                                workspace.settings.relativeGapTolerance,
                                                workspace.settings.gapsComputationMode)

    return min(threshold,workspace.settings.objectiveCutoff)
end
