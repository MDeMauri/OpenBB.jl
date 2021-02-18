# @Author: Massimo De Mauri <massimo>
# @Date:   2019-10-17T13:54:36+02:00
# @Email:  massimo.demauri@gmail.com
# @Filename: update_node.jl
# @Last modified by:   massimo
# @Last modified time: 2021-01-26T13:42:57+01:00
# @License: LGPL-3.0
# @Copyright: {{copyright}}

# this function updates a node to the current version
function update!(node::BBnode,workspace::BBworkspace)::Tuple{Bool,Bool}

	# check if there is anything to do at all
	if node.version == workspace.updatesRegister.masterVersion
		return (false,true)
	end


	# check if the register and the node version are coherent
	if node.version > workspace.updatesRegister.masterVersion
		@info (workspace.updatesRegister.masterVersion,workspace.updatesRegister.updates,node.version)
	end
	@assert node.version <= workspace.updatesRegister.masterVersion

	# use a boolean array to keep track of the updated variables
	updatedVars = falses(get_numVariables(workspace))
	updatedCnss = falses(get_numConstraints(workspace))

	# apply the all the updates necessary
	for updateID in node.version+1:workspace.updatesRegister.masterVersion
		feasible = workspace.updatesRegister.updates[updateID](node,workspace.updatesRegister.arguments[updateID]...,workspace,updatedVars,updatedCnss)
		if !feasible # the update can fail if the problem is discovered infeasible
			return (true,false)
		end
	end

	# declare the verision of the node to be the same of the master
	node.version = workspace.updatesRegister.masterVersion

	# finish collecting variables to preprocess
	for k in findall(updatedCnss)
		updatedVars[get_dependency(workspace.problem.cnsSet,k)] .= true
	end

	# preprocess
	return (true,preprocess_node!(node, workspace, findall(updatedVars),withBoundsPropagation=workspace.settings.withBoundsPropagation))
end


# dummy update
function do_nothing!(node::BBnode,workspace::BBworkspace,
					 newUpdatedVars::BitArray{1},newUpdatedCnss::BitArray{1})::Bool
	return true
end

# adapt the node to the insertion of constraints
function insert_constraints!(node::BBnode,position::Int64,newBounds::Tuple{Vector{Float},Vector{Float}},workspace::BBworkspace,
							 newUpdatedVars::BitArray{1},newUpdatedCnss::BitArray{1})::Bool
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
function remove_constraints!(node::BBnode,indices::Union{Vector{Int},UnitRange{Int}},workspace::BBworkspace,
						     newUpdatedVars::BitArray{1},newUpdatedCnss::BitArray{1})::Bool

	# delete the duals of the constraints to remove
	deleteat!(node.cnsDual,indices)
	deleteat!(node.cnsLoBs,indices)
	deleteat!(node.cnsUpBs,indices)
	#remove the constraints from newUpdatedCnss
	deleteat!(newUpdatedCnss,indices)
	return true
end

# changes the constraints order
function permute_constraints!(node::BBnode,permutation::Vector{Int64},workspace::BBworkspace,
							  newUpdatedVars::BitArray{1},newUpdatedCnss::BitArray{1})::Bool
	# check input
	@assert length(permutation) == length(node.cnsDual)
	@assert length(permutation) == length(node.cnsLoBs)
	@assert length(permutation) == length(node.cnsUpBs)
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
function update_bounds!(node::BBnode,varLoBs::Vector{Float},varUpBs::Vector{Float},
						cnsLoBs::Vector{Float},cnsUpBs::Vector{Float},workspace::BBworkspace,
						newUpdatedVars::BitArray{1},newUpdatedCnss::BitArray{1})::Bool
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
function insert_variables!(node::BBnode,position::Int64,newPrimal::Vector{Float},newBounds::Tuple{Vector{Float},Vector{Float}},workspace::BBworkspace,
						   newUpdatedVars::BitArray{1},newUpdatedCnss::BitArray{1})::Bool
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

function append_variables!(node::BBnode,newPrimal::Vector{Float},newBounds::Tuple{Vector{Float},Vector{Float}},workspace::BBworkspace,
						   newUpdatedVars::BitArray{1},newUpdatedCnss::BitArray{1})::Bool
	return insert_variables!(node,length(node.primal)+1,newPrimal,newBounds,
							newUpdatedVars,newUpdatedCnss)
end


# rounds the variable bounds (usually because they were marked as integral)
function round_variable_bounds!(node::BBnode,indices::Vector{Int64},workspace::BBworkspace,
	newUpdatedVars::BitArray{1},newUpdatedCnss::BitArray{1})::Bool
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


# remove the indicated variables
function remove_variables!(node::BBnode,indices::Vector{Int},workspace::BBworkspace,
					       newUpdatedVars::BitArray{1},newUpdatedCnss::BitArray{1})::Bool

	# check input
	@assert length(indices) == length(values)
	# collect info
	primalTolerance = workspace.settings.primalTolerance

	# remove the requested variables
	deleteat!(node.primal,indices)
	deleteat!(node.varLoBs,indices)
	deleteat!(node.varUpBs,indices)
	deleteat!(node.bndDual,indices)

	# fix the variables in the cuts
	remove_variables!(node.cutSet,indices,values)

	# apply the change also to the updated variables array
	deleteat!(newUpdatedVars,indices)

	return true
end


# fix value of the given variables to the given assignment
function fix_variables!(node::BBnode,indices::Vector{Int},values::Vector{Float},workspace::BBworkspace,
					    newUpdatedVars::BitArray{1},newUpdatedCnss::BitArray{1};
						removeFixedVariables::Bool=false)::Bool

	# check input
	@assert length(indices) == length(values)

	# collect info
	primalTolerance = workspace.settings.primalTolerance

	# check if the new assignment is feasible and fix values
	for (k,i) in enumerate(indices)
		# check feasibility
		if !(node.varLoBs[ind]-primalTolerance<=values[k]<=node.varLoBs[ind]+primalTolerance)
			return false
		end
	end

	# compute the change in objective
	# actually computing the lower bound may be expensive as it requires the
	# storage of the constraints and the objective, anyhow if the new values
	# are feasible, fixing the new value cannot improve the lower bound
 	node.objUpB = Inf

	# remove the useless variables
	if removeFixedVariables
		deleteat!(node.primal,indices)
		deleteat!(node.varLoBs,indices)
		deleteat!(node.varUpBs,indices)
		deleteat!(node.bndDual,indices)
		deleteat!(newUpdatedVars,indices)
	else
		node.primal[indices] = values
		node.varLoBs[indices] = values
		node.varUpBs[indices] = values
		newUpdatedVars[indices] .= true
	end

	# fix the variables in the cuts
	fix_variables!(node.cutSet,indices,values,removeFixedVariables=removeFixedVariables)


	return true
end



# Recomputation of objective bounds
function update_objLoB!(node::BBnode,problem::Problem,dualTolerance::Float)::Nothing
    # check the input
    @assert problem.objFun isa LinearObjective || problem.objFun isa QuadraticObjective || problem.objFun isa ConvexObjective
    @assert problem.cnsSet isa LinearConstraintSet || problem.cnsSet isa ConvexConstraintSet

	# detect cuts
    withCuts = nnz(node.cutSet.A)>0

    # compute the value of the dual constraints
    objGrd = evaluate_gradient(problem.objFun,node.primal)
    cnsJbc = evaluate_jacobian(problem.cnsSet,node.primal)
    dualCnssVal = objGrd + cnsJbc'*node.cnsDual + node.bndDual
    if withCuts
        dualCnssVal += node.cutSet.A'*node.cutDual
    end

    # check dual feasibility
    if maximum(abs.(dualCnssVal)) < dualTolerance
        # compute the value of the objective and the constraints
        objVal = evaluate(problem.objFun,node.primal)
        cnsVal = evaluate(problem.cnsSet,node.primal)

        # compute the dual objective
        node.objLoB = objVal + cnsVal'*node.cnsDual + node.primal'*node.bndDual +
                      -sum(weakInfMult.(min.(node.bndDual,0.0),node.varLoBs,dualTolerance*1e-2)) +
            	      -sum(weakInfMult.(max.(node.bndDual,0.0),node.varUpBs,dualTolerance*1e-2)) +
            	      -sum(weakInfMult.(min.(node.cnsDual,0.0),node.cnsLoBs,dualTolerance*1e-2)) +
            	      -sum(weakInfMult.(max.(node.cnsDual,0.0),node.cnsUpBs,dualTolerance*1e-2))

        if withCuts
            node.objLoB += -sum(weakInfMult.(min.(node.cutDual,0.0),node.cutSet.loBs,dualTolerance*1e-2)) +
        		           -sum(weakInfMult.(max.(node.cutDual,0.0),node.cutSet.upBs,dualTolerance*1e-2))
        end

        # check for numerical errors
        if isnan(node.objLoB)
            node.objLoB = -Inf
        end
    else
        node.objLoB = -Inf
    end
    return
end

# specialization (not necessary but more performant)
function update_objLoB!(node::BBnode,problem::Problem{<:LinearObjective,<:LinearConstraintSet},dualTolerance::Float)::Nothing

    # detect cuts
    withCuts = nnz(node.cutSet.A)>0

    # check dual feasibility
    dualCnssVal = problem.objFun.L + problem.cnsSet.A'*node.cnsDual + node.bndDual
    if withCuts
        dualCnssVal += node.cutSet.A'*node.cutDual
    end
    if maximum(abs.(dualCnssVal)) < dualTolerance
        # recompute objective bound
        node.objLoB = problem.objFun.c +
					  -sum(weakInfMult.(min.(node.bndDual,0.0),node.varLoBs,dualTolerance*1e-2)) +
                      -sum(weakInfMult.(max.(node.bndDual,0.0),node.varUpBs,dualTolerance*1e-2)) +
                      -sum(weakInfMult.(min.(node.cnsDual,0.0),node.cnsLoBs,dualTolerance*1e-2)) +
                      -sum(weakInfMult.(max.(node.cnsDual,0.0),node.cnsUpBs,dualTolerance*1e-2))
        if withCuts
            node.objLoB += -sum(weakInfMult.(min.(node.cutDual,0.0),node.cutSet.loBs,dualTolerance*1e-2)) +
                           -sum(weakInfMult.(max.(node.cutDual,0.0),node.cutSet.upBs,dualTolerance*1e-2))
        end

        # check for numerical errors
        if isnan(node.objLoB)
            node.objLoB = -Inf
        end
    else
        node.objLoB = -Inf
    end
    return
end

# specialization (not necessary but more performant and more robust)
function update_objLoB!(node::BBnode,problem::Problem{<:QuadraticObjective,<:LinearConstraintSet},dualTolerance::Float)::Nothing

    # detect cuts
    withCuts = nnz(node.cutSet.A)>0

    # check dual feasibility
    dualCnssVal = problem.objFun.Q*node.primal + problem.objFun.L + problem.cnsSet.A'*node.cnsDual + node.bndDual
    if withCuts
        dualCnssVal += node.cutSet.A'*node.cutDual
    end
    # try to use the quadratic term to remove the possible dual infeasibility
    if maximum(abs.(dualCnssVal)) > dualTolerance && !isnothing(problem.objFun.pInvQ)
        deltaPrimal = -problem.objFun.pInvQ*dualCnssVal
        dualCnssVal += problem.objFun.Q*deltaPrimal
        surrogatePrimal = node.primal + deltaPrimal
    else
        surrogatePrimal = node.primal
    end

    if maximum(abs.(dualCnssVal)) <= dualTolerance
        # recompute objective bound
        node.objLoB = problem.objFun.c +
					  -sum(weakInfMult.(min.(node.bndDual,0.0),node.varLoBs,dualTolerance*1e-2)) +
        	          -sum(weakInfMult.(max.(node.bndDual,0.0),node.varUpBs,dualTolerance*1e-2)) +
        	          -sum(weakInfMult.(min.(node.cnsDual,0.0),node.cnsLoBs,dualTolerance*1e-2)) +
        	          -sum(weakInfMult.(max.(node.cnsDual,0.0),node.cnsUpBs,dualTolerance*1e-2)) +
                      -0.5*surrogatePrimal'*problem.objFun.Q*surrogatePrimal
        if withCuts
            node.objLoB += -sum(weakInfMult.(min.(node.cutDual,0.0),node.cutSet.loBs,dualTolerance*1e-2)) +
                		   -sum(weakInfMult.(max.(node.cutDual,0.0),node.cutSet.upBs,dualTolerance*1e-2))
        end

        # check for numerical errors
        if isnan(node.objLoB)
            node.objLoB = -Inf
        end
    else
        node.objLoB = -Inf
    end
    return
end



function update_objUpB!(node::BBnode,problem::Problem,primalTolerance::Float)::Nothing

    # check primal feasibility
    cnsVal  = evaluate(problem.cnsSet,node.primal)
    cnsLoBs = node.cnsLoBs
    cnsUpBs = node.cnsUpBs
    if nnz(node.cutSet.A) > 0
        cnsVal  = vcat(cnsVal,evaluate(node.cutSet,node.primal))
        cnsLoBs = vcat(cnsLoBs,node.cutSet.loBs)
        cnsUpBs = vcat(cnsUpBs,node.cutSet.upBs)
    end

    if all(@. cnsLoBs - primalTolerance <= cnsVal <= cnsUpBs + primalTolerance)
        node.objUpB = evaluate(problem.objFun,node.primal)
	else
		node.objUpB = Inf
    end
   return
end

function update_objBounds!(node::BBnode,problem::Problem,primalTolerance::Float,dualTolerance::Float)::Nothing
	update_objLoB!(node,problem,dualTolerance), update_objUpB!(node,problem,primalTolerance)
	return
end
