# @Author: Massimo De Mauri <massimo>
# @Date:   2019-10-17T13:54:36+02:00
# @Email:  massimo.demauri@gmail.com
# @Filename: update_node.jl
# @Last modified by:   massimo
# @Last modified time: 2019-11-29T15:26:02+01:00
# @License: LGPL-3.0
# @Copyright: {{copyright}}

# this function updates a node to the current version
function update!(node::BBnode,workspace::BBworkspace{T1,T2,T3})::Bool where T1<:Problem where T2<:AbstractWorkspace where T3<:AbstractSharedMemory

	# check that the register is correct
	@assert node.version <= workspace.updatesRegister.masterVersion

	# use a boolean array to keep track of the updated variables
	updatedVars = falses(get_numVariables(workspace))
	updatedCnss = falses(get_numConstraints(workspace))

	# apply the all the updates necessary
	for updateID in node.version+1:workspace.updatesRegister.masterVersion
		feasible = workspace.updatesRegister.updates[updateID](node,workspace.updatesRegister.arguments[updateID]...,workspace,updatedVars,updatedCnss)
		if !feasible # the update can fail if the problem is discovered infeasible
			return false
		end
	end

	# declare the verision of the node to be the same of the master
	node.version = workspace.updatesRegister.masterVersion

	# finish collecting variables to preprocess
	for k in findall(updatedCnss)
		updatedVars[get_sparsity(workspace.problem.cnsSet,k)] .= true
	end
	# preprocess
	return preprocess!(node, workspace, findall(updatedVars))
end


# adapt the node to the insertion of constraints
function insert_constraints!(node::BBnode,position::Int64,newBounds::Tuple{Array{Float64,1},Array{Float64,1}},workspace::BBworkspace{T1,T2,T3},
							 newUpdatedVars::BitArray{1},newUpdatedCnss::BitArray{1})::Bool where T1<:Problem where T2<:AbstractWorkspace where T3<:AbstractSharedMemory
	# check input
	@assert length(newBounds[1]) == length(newBounds[2])
	numNewCnss = length(newBounds[1])
	# extend bounds
	splice!(node.cnsDual,position:position-1,zeros(numNewCnss))
	splice!(node.cnsLoBs,position:position-1,copy(newBounds[1]))
	# extend dual solution
	splice!(node.cnsUpBs,position:position-1,copy(newBounds[2]))
	# shift accordingly the newUpdatedCnss set (extending it if necessary)
	if position+numNewCnss > length(newUpdatedCnss)
		append!(newUpdatedCnss,falses(position+numNewCnss-length(newUpdatedCnss)))
	end
	@. newUpdatedCnss[position+numNewCnss:end] = newUpdatedCnss[position:end-numNewCnss]
	@. newUpdatedCnss[position:position+numNewCnss-1] = true

	return true
end

# adapt the node to the removal of constraints
function remove_constraints!(node::BBnode,indices::Array{Int64,1},workspace::BBworkspace{T1,T2,T3},
							 newUpdatedVars::BitArray{1},newUpdatedCnss::BitArray{1})::Bool where T1<:Problem where T2<:AbstractWorkspace where T3<:AbstractSharedMemory
	deleteat!(node.cnsDual,indices)
	deleteat!(node.cnsLoBs,indices)
	deleteat!(node.cnsUpBs,indices)
	#remove the constraints from newUpdatedCnss
	deleteat!(newUpdatedCnss,indices[findall(x->x<=length(newUpdatedCnss),indices)])
	return true
end

# changes the constraints order
function permute_constraints!(node::BBnode,permutation::Array{Int64,1},workspace::BBworkspace{T1,T2,T3},
							  newUpdatedVars::BitArray{1},newUpdatedCnss::BitArray{1})::Bool where T1<:Problem where T2<:AbstractWorkspace where T3<:AbstractSharedMemory
	# check input
	@assert length(permutation) == length(node.cnsDual)
	# permute the constraints bounds and dual solution
	permute!(node.cnsDual,permutation)
	permute!(node.cnsLoBs,permutation)
	permute!(node.cnsUpBs,permutation)
	# permute newUpdatedCnss (extending it if necessary)
	if length(permutation) > length(newUpdatedCnss)
		append!(newUpdatedCnss,falses(length(permutation)-length(newUpdatedCnss)))
	end
	newUpdatedCnss[1:length(permutation)] .= (newUpdatedCnss[1:length(permutation)])[permutation]
	return true
end


# changes variable and constraint bounds in the node
function update_bounds!(node::BBnode,varLoBs::Array{Float64,1},varUpBs::Array{Float64,1},
						cnsLoBs::Array{Float64,1},cnsUpBs::Array{Float64,1},workspace::BBworkspace{T1,T2,T3},
						newUpdatedVars::BitArray{1},newUpdatedCnss::BitArray{1})::Bool where T1<:Problem where T2<:AbstractWorkspace where T3<:AbstractSharedMemory
	# check input
	@assert length(varLoBs) == 0 || length(varLoBs) == length(node.varLoBs)
	@assert length(varUpBs) == 0 || length(varUpBs) == length(node.varUpBs)
	@assert length(cnsLoBs) == 0 || length(cnsLoBs) == length(node.cnsLoBs)
	@assert length(cnsUpBs) == 0 || length(cnsUpBs) == length(node.cnsUpBs)

	# extend newUpdatedCnss if necessary
	if length(cnsLoBs) > length(newUpdatedCnss)
		append!(newUpdatedCnss,falses(length(cnsLoBs)-length(newUpdatedCnss)))
	end

	# update the variable bounds
	if length(varLoBs) > 0
		for k in 1:length(node.varLoBs)
			if varLoBs[k] > node.varUpBs[k]
				return false
			elseif varLoBs[k] > node.varLoBs[k]
				node.varLoBs[k] = varLoBs[k]
				newUpdatedVars[k] = true
			end
		end
	end
	if length(varUpBs) > 0
		for k in 1:length(node.varUpBs)
			if varUpBs[k] < node.varLoBs[k]
				return false
			elseif varUpBs[k] < node.varUpBs[k]
				node.varUpBs[k] = varUpBs[k]
				newUpdatedVars[k] = true
			end
		end
	end

	# update the constraints bounds
	if length(cnsLoBs) > 0
		for k in 1:length(node.cnsLoBs)
			if cnsLoBs[k] > node.cnsUpBs[k]
				return false
			elseif cnsLoBs[k] > node.cnsLoBs[k]
				node.cnsLoBs[k] = cnsLoBs[k]
				newUpdatedCnss[k] = true
			end
		end
	end
	if length(cnsUpBs) > 0
		for k in 1:length(node.cnsUpBs)
			if cnsUpBs[k] < node.cnsLoBs[k]
				return false
			elseif cnsUpBs[k] < node.cnsUpBs[k]
				node.cnsUpBs[k] = cnsUpBs[k]
				newUpdatedCnss[k] = true
			end
		end
	end

	return true
end


# insert a new set of variables in the node data
function insert_variables!(node::BBnode,position::Int64,newPrimal::Array{Float64,1},newBounds::Tuple{Array{Float64,1},Array{Float64,1}},workspace::BBworkspace{T1,T2,T3},
						   newUpdatedVars::BitArray{1},newUpdatedCnss::BitArray{1})::Bool where T1<:Problem where T2<:AbstractWorkspace where T3<:AbstractSharedMemory
	# check input
	@assert length(newBounds[1]) == length(newBounds[2]) == length(newPrimal)
	numNewVars = length(newPrimal)
	# extend bounds arrays
	splice!(node.varLoBs,position:position-1,copy(newBounds[1]))
	splice!(node.varUpBs,position:position-1,copy(newBounds[2]))
	# extend primal/dual solutions
	splice!(node.primal,position:position-1,copy(newPrimal))
	splice!(node.bndDual,position:position-1,zeros(numNewVars))
	# shift accordingly the newUpdatedVars set
	@. newUpdatedVars[position+numNewVars:end] = newUpdatedVars[position:end-numNewVars]
	@. newUpdatedVars[position:position+numNewVars-1] = true

	return true
end


# rounds the variable bounds (usually because they were marked as integral)
function round_variable_bounds!(node::BBnode,indices::Array{Int64,1},workspace::BBworkspace{T1,T2,T3},
	newUpdatedVars::BitArray{1},newUpdatedCnss::BitArray{1})::Bool where T1<:Problem where T2<:AbstractWorkspace where T3<:AbstractSharedMemory
	# round bounds
	@. node.varLoBs[indices] = ceil(node.varLoBs[indices]-workspace.settings.primalTolerance)
	@. node.varUpBs[indices] = floor(node.varUpBs[indices]+workspace.settings.primalTolerance)
	# center primal solution
	@. node.primal = min(max(node.primal,node.varLoBs),node.varUpBs)
	# check the feasibility of new bounds
	for k in 1:length(node.varLoBs)
		if node.varLoBs[k] > node.varUpBs[k] + 2*workspace.settings.primalTolerance
			return false
		end
	end

	# denote updated variables
	newUpdatedVars[indices] .= true

	return true
end


# fix value of the given variables to the given assignment
function fix_variables!(node::BBnode,indices::Array{Int,1},values::Array{Float64,1},workspace::BBworkspace{T1,T2,T3},
					    newUpdatedVars::BitArray{1},newUpdatedCnss::BitArray{1})::Bool where T1<:Problem where T2<:AbstractWorkspace where T3<:AbstractSharedMemory

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

	# denote updated variables
	newUpdatedVars[indices] .= true

	return true
end
