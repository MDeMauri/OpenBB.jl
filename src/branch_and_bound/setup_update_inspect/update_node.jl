# @Author: Massimo De Mauri <massimo>
# @Date:   2019-10-17T13:54:36+02:00
# @Email:  massimo.demauri@gmail.com
# @Filename: update_node.jl
# @Last modified by:   massimo
# @Last modified time: 2019-11-18T11:25:02+01:00
# @License: LGPL-3.0
# @Copyright: {{copyright}}

# this function updates a node to the current version
function update!(node::BBnode,register::BBupdatesRegister)::Bool

    # check that the register is correct
    @assert node.version <= register.masterVersion

    # apply the all the updates necessary
    for updateID in node.version+1:register.masterVersion
        if register.updates[updateID](node,register.arguments[updateID]...) == false # the update can fail if the problem is discovered infeasible
			return false
		end
    end

	# declare the verision of the node to be the same of the master
	node.version = register.masterVersion

    return true
end


# adapt the node to the insertion of constraints
function insert_constraints!(node::BBnode,position::Int64,newBounds::Tuple{Array{Float64,1},Array{Float64,1}})::Bool
	@assert length(newBounds[1]) == length(newBounds[2])
	splice!(node.cnsDual,position:position-1,zeros(length(newBounds[1])))
    splice!(node.cnsLoBs,position:position-1,copy(newBounds[1]))
	splice!(node.cnsUpBs,position:position-1,copy(newBounds[2]))
	return true
end

# adapt the node to the removal of constraints
function remove_constraints!(node::BBnode,indices::Array{Int64,1})::Bool
    deleteat!(node.cnsDual,indices)
	deleteat!(node.cnsLoBs,indices)
	deleteat!(node.cnsUpBs,indices)
	return true
end

# changes the constraints order
function permute_constraints!(node::BBnode,permutation::Array{Int64,1})::Bool
	permute!(node.cnsDual,permutation)
	permute!(node.cnsLoBs,permutation)
	permute!(node.cnsUpBs,permutation)
	return true
end


# changes variable and constraint bounds in the node
function update_bounds!(node::BBnode,varLoBs::Array{Float64,1},varUpBs::Array{Float64,1},cnsLoBs::Array{Float64,1},cnsUpBs::Array{Float64,1})::Bool
	if length(varLoBs)>0 @. node.varLoBs = max(node.varLoBs,varLoBs) end
	if length(varUpBs)>0 @. node.varUpBs = min(node.varUpBs,varUpBs) end
	if length(cnsLoBs)>0 @. node.cnsLoBs = max(node.cnsLoBs,cnsLoBs) end
	if length(cnsUpBs)>0 @. node.cnsUpBs = min(node.cnsUpBs,cnsUpBs) end
	return true
end


# insert a new set of variables in the node data
function insert_variables!(node::BBnode,position::Int64,newPrimal::Array{Float64,1},newBounds::Tuple{Array{Float64,1},Array{Float64,1}})::Bool
	@assert length(newBounds[1]) == length(newBounds[2]) == length(newPrimal)
    splice!(node.varLoBs,position:position-1,copy(newBounds[1]))
    splice!(node.varUpBs,position:position-1,copy(newBounds[2]))
	splice!(node.primal,position:position-1,copy(newPrimal))
	splice!(node.bndDual,position:position-1,zeros(length(newPrimal)))
	return true
end


# rounds the variable bounds (usually because they were marked as integral)
function round_variable_bounds!(node::BBnode,indices::Array{Int64,1},tolerance::Float64)::Bool
	@. node.varLoBs[indices] =  ceil(node.varLoBs[indices]-tolerance)
	@. node.varUpBs[indices] = floor(node.varUpBs[indices]+tolerance)
	return true
end


# fix value of the given variables to the given assignment
function fix_variables!(node::BBnode,indices::Array{Int,1},values::Array{Float64,1},tolerance::Float64)::Bool

	# check if the assignment is feasible for the node
	for i in indices
		if !(node.varLoBs[i]-tolerance < values[i] < node.varUpBs[i] + tolerance)
			return false
		end
	end
	# fix values
	@. node.varLoBs[indices] = max(node.varLoBs[indices],values)
	@. node.varUpBs[indices] = min(node.varUpBs[indices],values)
	@. node.primal[indices] = values
	return true
end
