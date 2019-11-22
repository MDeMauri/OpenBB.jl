# @Author: Massimo De Mauri <massimo>
# @Date:   2019-10-17T13:54:36+02:00
# @Email:  massimo.demauri@gmail.com
# @Filename: update_node.jl
# @Last modified by:   massimo
# @Last modified time: 2019-11-22T10:35:58+01:00
# @License: LGPL-3.0
# @Copyright: {{copyright}}

# this function updates a node to the current version
function update!(node::BBnode,workspace::BBworkspace{T1,T2,T3})::Bool where T1<:Problem where T2<:AbstractWorkspace where T3<:AbstractSharedMemory

    # check that the register is correct
    @assert node.version <= workspace.updatesRegister.masterVersion

    # apply the all the updates necessary
    for updateID in node.version+1:workspace.updatesRegister.masterVersion
        if workspace.updatesRegister.updates[updateID](node,workspace.updatesRegister.arguments[updateID]...,workspace) == false # the update can fail if the problem is discovered infeasible
			return false
		end
    end

	# declare the verision of the node to be the same of the master
	node.version = workspace.updatesRegister.masterVersion

    return true
end


# adapt the node to the insertion of constraints
function insert_constraints!(node::BBnode,position::Int64,newBounds::Tuple{Array{Float64,1},Array{Float64,1}},
							 workspace::BBworkspace{T1,T2,T3})::Bool where T1<:Problem where T2<:AbstractWorkspace where T3<:AbstractSharedMemory
	# check input
	@assert length(newBounds[1]) == length(newBounds[2])
	# extend bounds
	splice!(node.cnsDual,position:position-1,zeros(length(newBounds[1])))
    splice!(node.cnsLoBs,position:position-1,copy(newBounds[1]))
	# extend dual solution
	splice!(node.cnsUpBs,position:position-1,copy(newBounds[2]))
	#TODO preprocessing
	return true
end

# adapt the node to the removal of constraints
function remove_constraints!(node::BBnode,indices::Array{Int64,1},
							 workspace::BBworkspace{T1,T2,T3})::Bool where T1<:Problem where T2<:AbstractWorkspace where T3<:AbstractSharedMemory
    deleteat!(node.cnsDual,indices)
	deleteat!(node.cnsLoBs,indices)
	deleteat!(node.cnsUpBs,indices)
	return true
end

# changes the constraints order
function permute_constraints!(node::BBnode,permutation::Array{Int64,1},
							  workspace::BBworkspace{T1,T2,T3})::Bool where T1<:Problem where T2<:AbstractWorkspace where T3<:AbstractSharedMemory
	# check input
	@assert length(permutation) == length(node.cnsDual)
	# permute the constraints bounds and dual solution
	permute!(node.cnsDual,permutation)
	permute!(node.cnsLoBs,permutation)
	permute!(node.cnsUpBs,permutation)
	return true
end


# changes variable and constraint bounds in the node
function update_bounds!(node::BBnode,varLoBs::Array{Float64,1},varUpBs::Array{Float64,1},
						cnsLoBs::Array{Float64,1},cnsUpBs::Array{Float64,1},
						workspace::BBworkspace{T1,T2,T3})::Bool where T1<:Problem where T2<:AbstractWorkspace where T3<:AbstractSharedMemory
	# check input
	@assert length(varLoBs) == 0 || length(varLoBs) == length(node.varLoBs)
	@assert length(varUpBs) == 0 || length(varUpBs) == length(node.varUpBs)
	@assert length(cnsLoBs) == 0 || length(cnsLoBs) == length(node.cnsLoBs)
	@assert length(cnsUpBs) == 0 || length(cnsUpBs) == length(node.cnsUpBs)
	# update the bounds
	if length(varLoBs)>0 @. node.varLoBs = max(node.varLoBs,varLoBs) end
	if length(varUpBs)>0 @. node.varUpBs = min(node.varUpBs,varUpBs) end
	if length(cnsLoBs)>0 @. node.cnsLoBs = max(node.cnsLoBs,cnsLoBs) end
	if length(cnsUpBs)>0 @. node.cnsUpBs = min(node.cnsUpBs,cnsUpBs) end
	# check the feasibility of new bounds
	for k in 1:length(node.varLoBs)
		if node.varLoBs[k] > node.varUpBs[k] + 2*workspace.settings.primalTolerance
			return false
		end
	end
	for k in 1:length(node.cnsLoBs)
		if node.cnsLoBs[k] > node.cnsUpBs[k] + 2*workspace.settings.primalTolerance
			return false
		end
	end
	#TODO preprocessing
	return true
end


# insert a new set of variables in the node data
function insert_variables!(node::BBnode,position::Int64,newPrimal::Array{Float64,1},newBounds::Tuple{Array{Float64,1},Array{Float64,1}},
						   workspace::BBworkspace{T1,T2,T3})::Bool where T1<:Problem where T2<:AbstractWorkspace where T3<:AbstractSharedMemory
	# check input
	@assert length(newBounds[1]) == length(newBounds[2]) == length(newPrimal)
	# extend bounds arrays
    splice!(node.varLoBs,position:position-1,copy(newBounds[1]))
    splice!(node.varUpBs,position:position-1,copy(newBounds[2]))
	# extend primal/dual solutions
	splice!(node.primal,position:position-1,copy(newPrimal))
	splice!(node.bndDual,position:position-1,zeros(length(newPrimal)))
	return true
end


# rounds the variable bounds (usually because they were marked as integral)
function round_variable_bounds!(node::BBnode,indices::Array{Int64,1},
								workspace::BBworkspace{T1,T2,T3})::Bool where T1<:Problem where T2<:AbstractWorkspace where T3<:AbstractSharedMemory
	# round bounds
	@. node.varLoBs[indices] =  ceil(node.varLoBs[indices]-workspace.settings.primalTolerance)
	@. node.varUpBs[indices] = floor(node.varUpBs[indices]+workspace.settings.primalTolerance)
	# center primal solution
	@. node.primal = min(max(node.primal,node.varLoBs),node.varUpBs)
	# check the feasibility of new bounds
	for k in 1:length(node.varLoBs)
		if node.varLoBs[k] > node.varUpBs[k] + 2*workspace.settings.primalTolerance
			return false
		end
	end
	#TODO preprocessing
	return true
end


# fix value of the given variables to the given assignment
function fix_variables!(node::BBnode,indices::Array{Int,1},values::Array{Float64,1},
						workspace::BBworkspace{T1,T2,T3})::Bool where T1<:Problem where T2<:AbstractWorkspace where T3<:AbstractSharedMemory
	# check input
	@assert length(indices) == length(values)
	# check if the new assignment is feasible and fix values
	for (k,i) in enumerate(indices)
		node.varLoBs[i] = max(node.varLoBs[i],values[k])
		node.varUpBs[i] = min(node.varUpBs[i],values[k])
		node.primal[i] = values[k]
		if node.varLoBs[i] > node.varUpBs[i] + 2*workspace.settings.primalTolerance
			return false
		end
	end
	#TODO preprocessing
	return true
end
