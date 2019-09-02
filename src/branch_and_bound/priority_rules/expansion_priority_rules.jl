# @Author: Massimo De Mauri <massimodemauri>
# @Date:   2019-02-06T16:19:45+01:00
# @Email:  massimo.demauri@gmaileftNode.com
# @Project: OpenBB
# @Filename: nodes_priority_functions.jl
# @Last modified by:   massimo
# @Last modified time: 2019-09-02T15:07:58+02:00
# @License: LGPL-3.0
# @Copyright: {{copyright}}


# wrappers for nodes priority rules
function expansion_priority_rule(functionTuple::Tuple,leftNode::BBnode,rightNode::BBnode,status::BBstatus)::Bool
    @assert functionTuple[1] isa Function
    return functionTuple[1](leftNode,rightNode,status,functionTuple[2:end]...)
end


# priority functions for nodes
function lower_objective(leftNode::BBnode,rightNode::BBnode,status::BBstatus)::Bool
    if leftNode.objVal <= rightNode.objVal
        return true
    else
        return false
    end
end

function higher_objective(leftNode::BBnode,rightNode::BBnode,status::BBstatus)::Bool
    if leftNode.objVal >= rightNode.objVal
        return true
    else
        return false
    end
end

function lower_avgAbsFractionality(leftNode::BBnode,rightNode::BBnode,status::BBstatus)::Bool
    if leftNode.avgAbsFrac <= rightNode.avgAbsFrac
        return true
    else
        return false
    end
end

function higher_avgAbsFractionality(leftNode::BBnode,rightNode::BBnode,status::BBstatus)::Bool
    if leftNode.avgAbsFrac >= rightNode.avgAbsFrac
        return true
    else
        return false
    end
end




# priority functions for nodes
function lower_pseudoObjective(leftNode::BBnode,rightNode::BBnode,status::BBstatus)::Bool
    if leftNode.pseudoObjective <= rightNode.pseudoObjective
        return true
    else
        return false
    end
end

function higher_pseudoObjective(leftNode::BBnode,rightNode::BBnode,status::BBstatus)::Bool
    if leftNode.pseudoObjective >= rightNode.pseudoObjective
        return true
    else
        return false
    end
end
