# @Author: Massimo De Mauri <massimodemauri>
# @Date:   2019-02-06T16:19:45+01:00
# @Email:  massimo.demauri@gmail.com
# @Project: OpenBB
# @Filename: nodes_priority_functions.jl
# @Last modified by:   massimo
# @Last modified time: 2021-01-28T15:18:46+01:00
# @License: LGPL-3.0
# @Copyright: {{copyright}}


# wrappers for nodes priority rules
function expansion_priority(node::BBnode,workspace::BBworkspace)::Float
    return OpenBB.eval(workspace.settings.expansionPriorityRule[1])(node,workspace,workspace.settings.expansionPriorityRule[2:end]...)
end


# priority functions for nodes
function lower_objective(node::BBnode,workspace::BBworkspace)::Float
    return -node.objLoB
end

function higher_objective(node::BBnode,workspace::BBworkspace)::Float
    return node.objUpB
end

function lower_avgFractionality(node::BBnode,workspace::BBworkspace)::Float
    return -node.avgFractionality
end

function higher_avgFractionality(node::BBnode,workspace::BBworkspace)::Float
    return node.avgFractionality
end


function lower_pseudoObjective(node::BBnode,workspace::BBworkspace)::Float
    return -(node.objLoB+node.pseudoCost)
end

function higher_pseudoObjective(node::BBnode,workspace::BBworkspace)::Float
    return node.objLoB+node.pseudoCost
end


function lower_dualTradeoff(node::BBnode,workspace::BBworkspace,gainOfDual::Float=1.0)::Float
    if gainOfDual == 0.0
        return node.objLoB
    else
        bndDualCost = sum(weakInfMult.(min.(node.bndDual,0.0),node.varLoBs,workspace.settings.dualTolerance)) +
                      sum(weakInfMult.(max.(node.bndDual,0.0),node.varUpBs,workspace.settings.dualTolerance))
        return -node.objLoB - gainOfDual*bndDualCost
    end
end

function higher_dualTradeoff(node::BBnode,workspace::BBworkspace,gainOfDual::Float=1.0)::Float
    if gainOfDual == 0.0
        return -node.objLoB
    else
        bndDualCost = sum(weakInfMult.(min.(node.bndDual,0.0),node.varLoBs,workspace.settings.dualTolerance)) +
                      sum(weakInfMult.(max.(node.bndDual,0.0),node.varUpBs,workspace.settings.dualTolerance))
        return node.objLoB + gainOfDual*bndDualCost
    end
end
