# @Author: Massimo De Mauri <massimo>
# @Date:   2019-06-07T18:58:51+02:00
# @Email:  massimo.demauri@gmail.com
# @Filename: pseudo_costs_initialization.jl
# @Last modified by:   massimo
# @Last modified time: 2019-09-02T15:08:06+02:00
# @License: LGPL-3.0
# @Copyright: {{copyright}}

# wrappers for initialization functions
function initialize_pseudoCosts!(functionTuple::Tuple,pseudoCosts::Tuple{Array{Float64,2},Array{Int,2}},refNode::BBnode)::Nothing
    @assert functionTuple[1] isa Function
    return functionTuple[1](pseudoCosts,refNode,functionTuple[2:end]...)
end


# actual initialization functions
function initialize_to_constant!(pseudoCosts::Tuple{Array{Float64,2},Array{Int,2}},refNode::BBnode,value::Float64)::Nothing
    @. pseudoCosts[1] = value
    @. pseudoCosts[2] = 0
    return
end

function initialize_with_root_objective!(pseudoCosts::Tuple{Array{Float64,2},Array{Int,2}},refNode::BBnode,proportionalityConstant::Float64)::Nothing
    @. pseudoCosts[1] = max(proportionalityConstant*abs(refNode.objUpB),1e-4)
    @. pseudoCosts[2] = 0
    return
end

function initialize_with_strong_branching!(pseudoCosts::Tuple{Array{Float64,2},Array{Int,2}},refNode::BBnode)::Nothing
    #TODO
    @error "Not Implemented Yet"
    return
end
