# @Author: Massimo De Mauri <massimo>
# @Date:   2019-05-20T10:04:17+02:00
# @Email:  massimo.demauri@gmail.com
# @Filename: update_nodes.jl
# @Last modified by:   massimo
# @Last modified time: 2021-02-14T19:22:02+01:00
# @License: LGPL-3.0
# @Copyright: {{copyright}}

# it marks the workspace as outdated
function make_outdated!(workspace::BBworkspace;infeasiblesToRecover::Bool=false,invalidLowerBounds::Bool=false)::Nothing
	workspace.outdated = true
	make_outdated!(workspace.subsolverWS)
	make_outdated!(workspace.heuristicsWS)
	if infeasiblesToRecover
		workspace.infeasiblesToRecover = true
	end
	if invalidLowerBounds
		workspace.invalidLowerBounds = true
	end
	return
end

# ...
function reset_global_info!(workspace::BBworkspace)::Nothing

	# reset the global info
	workspace.sharedMemory.globalObjUpB[1] = (Inf,Int8(0))
	@. workspace.sharedMemory.localObjLoBs = -Inf
	@. workspace.sharedMemory.stats = 0
	@. workspace.sharedMemory.arrestable = false

	return
end


#
function update_sharedMemory!(workspace::BBworkspace)::Nothing

	# update the communtication Channels
	if workspace.sharedMemory isa BBsharedMemory
		numVars = get_numVariables(workspace)
		numCnss = get_numConstraints(workspace)
		numDscVars = get_numDiscrete(workspace)
		currentMaxSerialSize = maxSerialSize_BBnode(numVars,numCnss,workspace.settings.maxNumberOfLocalCuts + is_mixedBinary(workspace.problem.varSet))

		# construct new communication Channels if needed
		if currentMaxSerialSize > get_size(workspace.sharedMemory.inputChannel)
			communicationChannels = Vector{BBnodeChannel}(undef,workspace.settings.numProcesses)
			for k in 1:workspace.settings.numProcesses
				communicationChannels[k] = BBnodeChannel(currentMaxSerialSize)
			end
			workspace.sharedMemory.inputChannel = communicationChannels[1]
			workspace.sharedMemory.outputChannel = communicationChannels[2]
			@sync for k in 2:workspace.settings.numProcesses
				@async if k < workspace.settings.numProcesses
					remotecall_fetch(OpenBB.eval,k,:(workspace.sharedMemory.inputChannel = $(communicationChannels[k]);
												   workspace.sharedMemory.outputChannel = $(communicationChannels[k+1]);
												   nothing))
				else
					remotecall_fetch(OpenBB.eval,k,:(workspace.sharedMemory.inputChannel = $(communicationChannels[k]);
												   workspace.sharedMemory.outputChannel = $(communicationChannels[1]);
												   nothing))
				end
			end
		end
	end
	return
end


#
function update!(workspace::BBworkspace;updateAllNodes::Bool=false,subsolverIterations=0,localOnly::Bool=false)::Nothing

    @sync if !localOnly && !(workspace.sharedMemory isa NullSharedMemory)

        # call the local version of the function on the remote workers
        for p in workers()[1:workspace.settings.numProcesses-1]
            @async remotecall_fetch(OpenBB.eval,p,:(OpenBB.update!(workspace,subsolverIterations=$subsolverIterations,
																			 localOnly=true)))
        end

        # call the local version of the function on the current process
        update!(workspace,subsolverIterations=subsolverIterations,localOnly=true)
		# reset the information shared among the processes
		reset_global_info!(workspace)
		# adapt the shared memory to the new problem
		update_sharedMemory!(workspace)
    else

		# adapt the node pools to the changes
		if workspace.infeasiblesToRecover
			workspace.tree.active = vcat(workspace.tree.infeasibles,workspace.tree.suboptimals,workspace.tree.active,workspace.tree.solutions)
			empty!(workspace.tree.solutions)
			empty!(workspace.tree.infeasibles)
			empty!(workspace.tree.suboptimals)
		else
			workspace.tree.active = vcat(workspace.tree.suboptimals,workspace.tree.active,workspace.tree.solutions)
			empty!(workspace.tree.solutions)
			empty!(workspace.tree.suboptimals)
		end

		# recompute lower bounds for the nodes
		if workspace.invalidLowerBounds

			# if needed, create a temporary subsolver for the recomputation of the lower-bounds
			if subsolverIterations != 0
				_subsolverSettings = deepcopy(workspace.subsolverWS.settings)
				if subsolverIterations != -1
					set_iterationsLimit!(_subsolverSettings,subsolverIterations)
				end
				_subsolver = setup(workspace.problem,_subsolverSettings)
			end

			# recompute lower bound for nodes
			workspace.status.objLoB = Inf
			for (k,(node,~)) in enumerate(workspace.tree.active)

				# update node
				(feasible,~) = update!(node,workspace)

				# recompute bounds
				if !feasible
					node.objLoB = node.objUpB = Inf
				else
					if subsolverIterations == 0
						update_objLoB!(node,workspace.problem,10.0*workspace.settings.dualTolerance)
						node.objUpB = Inf
					else
						solve!(node,_subsolver)
					end
					# recompute node priority
					workspace.tree.active[k] = (node,expansion_priority(node,workspace))
					# update the general lower bound
					workspace.status.objLoB = min(workspace.status.objLoB,node.objLoB)
				end
			end

			# remove the infeasible nodes from the active pool
			toRemove = findall(nodeTuple->nodeTuple[1].objLoB==Inf,workspace.tree.active)
			if workspace.settings.conservativismLevel == 2
				append!(workspace.tree.infeasibles,workspace.tree.active[toRemove])
			end
			deleteat!(workspace.tree.active,toRemove)

		elseif updateAllNodes

			workspace.status.objLoB = Inf
			for (k,(node,~)) in enumerate(workspace.tree.active)

				# update node
				(feasible,~) = update!(node,workspace)

				# update node
				if !feasible
					node.objLoB = node.objUpB = Inf
				else
					# recompute node priority
					workspace.tree.active[k] = (node,expansion_priority(node,workspace))
					# update the general lower bound
					workspace.status.objLoB = min(workspace.status.objLoB,node.objLoB)
				end
			end

			# remove the infeasible nodes from the active pool
			toRemove = findall(nodeTuple->nodeTuple[1].objLoB==Inf,workspace.tree.active)
			if workspace.settings.conservativismLevel == 2
				append!(workspace.tree.infeasibles,workspace.tree.active[toRemove])
			end
			deleteat!(workspace.tree.active,toRemove)
		else
			workspace.status.objLoB = Inf
			for (node,~) in workspace.tree.active
				workspace.status.objLoB = min(workspace.status.objLoB,node.objLoB)
			end
		end

		# re-sort the active queue
		sort!(workspace.tree.active)

		# adapt the status to the changes
		workspace.status.objUpB = Inf
		workspace.status.absoluteGap = Inf
		workspace.status.relativeGap = Inf
		workspace.status.numSolutions = 0
		workspace.status.reliable = true
		workspace.status.cutoffActive = false
		workspace.status.blackListActive = false
		workspace.status.description = "interrupted"

		# mark the workspace as up to date
		workspace.outdated = false
		workspace.infeasiblesToRecover = false
		workspace.invalidLowerBounds = false
    end

    return
end


# eliminates all the generated nodes from the workspace
function clear!(workspace::BBworkspace;localOnly::Bool=false)::Nothing

    @sync if !localOnly && !(workspace.sharedMemory isa NullSharedMemory)

        # call function on the remote workers
        for p in workers()[1:workspace.settings.numProcesses-1]
            @async remotecall_fetch(OpenBB.eval,p,:(OpenBB.clear!(workspace,localOnly=true)))
        end

        # call function on the main process
        clear!(workspace,localOnly=true)
		# reset the information shared among the processes
		reset_global_info!(workspace)

    else
		# empty the tree
        empty!(workspace.tree)
        # reset the status
        workspace.status = BBstatus()
		workspace.status.objUpB = Inf
		workspace.status.objLoB = Inf
		workspace.status.description = "empty"
		workspace.status.cutoffActive = false
		workspace.status.blackListActive = false
    end

    return
end



# return the workspace to the initial state
function reset!(workspace::BBworkspace;localOnly::Bool=false)::Nothing

    @sync if !localOnly && !(workspace.sharedMemory isa NullSharedMemory)

        # remove all nodes in the remote workers
        for p in workers()[1:workspace.settings.numProcesses-1]
            @async remotecall_fetch(OpenBB.eval,p,:(OpenBB.clear!(workspace,localOnly=true)))
        end

        # call the local version of the function on the main process
        reset!(workspace,localOnly=true)
		# reset the information shared among the processes
		reset_global_info!(workspace)

    else
		# collect some data for the BBworkspace
	    numVars = get_numVariables(workspace)
		numCnss = get_numConstraints(workspace)
		numDscVars = get_numDiscrete(workspace)
        # eliminate all the generated nodes
        clear!(workspace,localOnly=true)
		# reset the updatesRegister
		reset!(workspace.updatesRegister)
		# set the status to default
		workspace.status = BBstatus()
		# re-insert the root node into the queue
		if myid() == 1
	        push_node!(workspace,:active,BBroot(workspace))
		end
    end

    return
end



# this function is used to eliminate nodes using a selection function
function Base.filter!(selector::Function,workspace::BBworkspace;suppressErrors::Bool=false,localOnly::Bool=false)::Nothing

	if myid()==1 && !suppressErrors && workspace.settings.conservativismLevel == 0
		@warn "In order to correctly filter away nodes in the tree \"conservativismLevel\" should be set to at least 1"
	end

	# store changes
	oldObjUpB = newObjLoB = Inf

	if !localOnly && !(workspace.sharedMemory isa NullSharedMemory)

		@sync if true # fake "if" for synchronization
			# call the function in the remote workers
	        for p in workers()[1:workspace.settings.numProcesses-1]
	            @async remotecall_fetch(OpenBB.eval,p,:(filter!($selector,workspace,suppressErrors=true,localOnly=true)))
	        end

			# call the local version of the function in the current process
			filter!(selector,workspace,suppressErrors=true,localOnly=true)
		end

		if myid() == 1

			# collect info
			numProcesses = workspace.settings.numProcesses
			remoteObjUpBs = Vector{Float}(undef,numProcesses); remoteObjUpBs[1] = workspace.status.objUpB
			remoteNumSols = Vector{Float}(undef,numProcesses); remoteNumSols[1] = workspace.status.numSolutions
			@sync for (k,id) in enumerate(workers()[1:numProcesses])
				@async (remoteObjUpBs[k+1],remoteNumSols[k+1]) =
					remotecall_fetch(OpenBB.eval,p,:((workspace.status.objUpB,workspace.status.numSolutions)))
			end

			# re-compute global info
			oldObjUpB = workspace.sharedMemory.globalObjUpB[1][1]
			minIdx = argmin(remoteObjUpBs)
			if remoteObjUpBs[minIdx] < Inf
				workspace.sharedMemory.globalObjUpB[1] = (remoteObjUpBs[minIdx],minIdx)
				workspace.status.objUpB = remoteObjUpBs[minIdx]
			else
				workspace.sharedMemory.globalObjUpB[1] = (Inf,0)
				workspace.status.objUpB = Inf
			end
			workspace.sharedMemory.stats[1] = workspace.status.numSolutions = sum(remoteNumSols)
			newObjLoB = workspace.sharedMemory.globalObjUpB[1][1]


			# recompute optimality gaps
			(workspace.status.absoluteGap,workspace.status.relativeGap) =
				compute_gaps(min(workspace.status.objLoB,objNodeInOutputChannel),workspace.status.objUpB,workspace.settings.gapsComputationMode)

		end

	else

		# remove the nodes from the tree
		filter!(selector,workspace.tree)

		# recompute upper bound and number of solutions
		oldObjUpB = workspace.status.objUpB
		if !isempty(workspace.tree.solutions)
			workspace.status.objUpB = workspace.tree.solutions[end][1].objUpB
			workspace.status.numSolutions = length(workspace.tree.solutions)
		else
			workspace.status.objUpB = Inf
			workspace.status.numSolutions = 0
		end
		newObjUpB = workspace.status.objUpB

		# recompute lower_bound
		newObjLoB = min(workspace.status.objUpB,workspace.settings.objectiveCutoff)
		for i in length(workspace.tree.active):-1:1
			# if we have unreliable problems in the tree.active we cannot update the lower bound
			if !workspace.tree.active[i][1].reliable
				newObjLoB = workspace.status.objLoB
				break
			elseif newObjLoB > workspace.tree.active[i][1].objLoB
				newObjLoB = workspace.tree.active[i][1].objLoB
			end
		end
		workspace.status.objLoB = newObjLoB

		# recompute optimality gaps
		(workspace.status.absoluteGap,workspace.status.relativeGap) =
			compute_gaps(workspace.status.objLoB,workspace.status.objUpB,workspace.settings.gapsComputationMode)

		# communicate the local changes
		if !(workspace.sharedMemory isa NullSharedMemory)
			workspace.sharedMemory.localObjLoBs[myid()] = workspace.status.objLoB
		end
	end

	# mark the workspaces as outdated if necessary
	if  myid()==1 && oldObjUpB != newObjUpB
		if !(workspace.sharedMemory isa NullSharedMemory)
			@sync for (k,id) in enumerate(workers()[1:numProcesses])
				@async remotecall_fetch(OpenBB.eval,p,:(make_outdated!(workspace)))
			end
		end
		make_outdated!(workspace,infeasiblesToRecover=false,invalidLowerBounds=false)
	end

	return
end


## Optimal Control Specific Functions

# this function inserts the given set of new constraints into the OpenBB workspace
# preserving the sparsity pattern typical of optimal control problems
function insert_constraints_oc!(workspace::BBworkspace,constraintSet::ConstraintSet;
                                suppressErrors::Bool=false,localOnly::Bool=false)::Nothing

	# preliminary check
	@assert all(workspace.settings.optimalControlInfo.>-1)

	# collect info
	numVarsPerStep = sum(workspace.settings.optimalControlInfo)

    # append the constraints in the workspace
    append_constraints!(workspace,constraintSet,suppressErrors=suppressErrors,localOnly=localOnly)

    # sort the constraints in the workspace
    sort_constraints_oc!(workspace,localOnly=localOnly)

    return
end



# this function sorts the constraints in the workspace in order to ensure the
# sparsity pattern typical of optimal control problems
function sort_constraints_oc!(workspace::BBworkspace;suppressErrors::Bool=false,localOnly::Bool=false)::Nothing

  # check input
  @assert all(workspace.settings.optimalControlInfo .>= 0)
  numVarsPerStep = sum(workspace.settings.optimalControlInfo)

  @sync if !localOnly && !(workspace.sharedMemory isa NullSharedMemory)
       # call the local version of the function on the remote workers
       for p in workers()[1:workspace.settings.numProcesses-1]
           @async remotecall_fetch(OpenBB.eval,p,:(
			 	OpenBB.sort_constraints_oc!(workspace,suppressErrors=$suppressErrors,localOnly=true)
			 ))
       end
       # call the local version of the function on the main process
       sort_constraints_oc!(workspace,suppressErrors=suppressErrors,localOnly=true)

  else
	    # sort constraints according to the shooting structure
		if !issorted_oc(workspace.problem.cnsSet,workspace.settings.optimalControlInfo)
			permutation = sortperm_oc(workspace.problem.cnsSet,workspace.settings.optimalControlInfo)
			permute_constraints!(workspace,permutation,localOnly=true)
		end

	    # mark the subsolver workspace as outdated
	    make_outdated!(workspace.subsolverWS)
	end

    return
end
