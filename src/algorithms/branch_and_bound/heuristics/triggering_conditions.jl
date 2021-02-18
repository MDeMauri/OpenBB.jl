# @Author: Massimo De Mauri <massimo>
# @Date:   2020-12-17T14:44:27+01:00
# @Email:  massimo.demauri@gmail.com
# @Filename: triggering_conditions.jl
# @Last modified by:   massimo
# @Last modified time: 2021-01-24T19:48:22+01:00
# @License: LGPL-3.0
# @Copyright: {{copyright}}



####################### Triggering Conditions #################################
# check satisfaction of triggering condition
function check_heuristicsTriggerCondition(node::BBnode,workspace::BBworkspace)::Bool
    return OpenBB.eval(workspace.settings.heuristicsTriggerCondition[1])(node,workspace,workspace.settings.heuristicsTriggerCondition[2:end]...)
end

# never trigger heuristics
function never_trigger(node::BBnode,workspace::BBworkspace)::Bool
    return false
end

# check if the node respects the triggering condition
function trigger_on_fractionality(node::BBnode,workspace::BBworkspace,fractionalityThreshold::Float)::Bool
    if node.avgFractionality < fractionalityThreshold
        return true
    else
        return false
   end
end

function trigger_on_optimality_and_fractionality(node::BBnode,workspace::BBworkspace,optimalityThreshold::Float,fractionalityThreshold::Float)::Bool
    absoluteGap,relativeGap = compute_gaps(node.objLoB,workspace.status.objUpB,workspace.settings.gapsComputationMode)
    if relativeGap < optimalityThreshold && node.avgFractionality < fractionalityThreshold
        return true
    else
        return false
   end
end

function trigger_on_optimality_or_fractionality(node::BBnode,workspace::BBworkspace,optimalityThreshold::Float,fractionalityThreshold::Float)::Bool
    absoluteGap,relativeGap = compute_gaps(node.objLoB,workspace.status.objUpB,workspace.settings.gapsComputationMode)
    if relativeGap < optimalityThreshold || node.avgFractionality < fractionalityThreshold
        return true
    else
        return false
   end
end


function trigger_on_pseudoObjective(node::BBnode,workspace::BBworkspace)::Bool
    expectedObjVal = node.objLoB+node.pseudoCost
    # use the expected objective value to predict if the rounded node may lead to a new incumbent
    if expectedObjVal > get_suboptimality_threshold(workspace) || workspace.status.objLoB == -Inf
        return false
    else
        absoluteGap,relativeGap = compute_gaps(workspace.status.objLoB,expectedObjVal,workspace.settings.gapsComputationMode)
        if absoluteGap <= workspace.settings.absoluteGapTolerance ||
           relativeGap <= workspace.settings.relativeGapTolerance
            return true
        else
            return false
       end
   end
end
