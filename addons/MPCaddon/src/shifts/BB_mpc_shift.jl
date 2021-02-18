# @Author: Massimo De Mauri <massimo>
# @Date:   2021-01-13T16:46:45+01:00
# @Email:  massimo.demauri@protonmail.com
# @Filename: mpc_shift.jl
# @Last modified by:   massimo
# @Last modified time: 2021-02-16T13:58:01+01:00
# @License: LGPL-3.0
# @Copyright: {{copyright}}
#

function mpc_shift!(workspace::BBworkspace,shiftSteps::Int,
					newTailCosts::ObjectiveFunction,newTailConstraints::ConstraintSet,newTailVariables::VariableSet=EmptyVariableSet();
					referenceSolution::Array{Float64,1}=Float64[],measuredState::Array{Float64,1}=Float64[],
					mode::Symbol=:fullRTI,approximateBellmanConstant::Float=Inf,subsolverIterations::Int=0,
					suppressErrors::Bool=false,localOnly::Bool=false)::Nothing

	if workspace.settings.verbose
		println("BB mpc_shift!, mode = ",mode)
	end

	# check workspace
	if all(workspace.settings.optimalControlInfo .== -1)
		raise(ArgumentError("mpc_shift!: The workspace has no optimal control info."))
	else
		@assert all(workspace.settings.optimalControlInfo .>= 0)
	end

	# check the problem
	if workspace.problem.objFun.c != 0
		error("HBB mpc_shift!: The algorithm cannot handle objective function with constant parts, please reformulate the problem (possibly using an algebraic variable)")
	end


	# collect info
	numStates = workspace.settings.optimalControlInfo[1]
	numAlgebraic = workspace.settings.optimalControlInfo[2]
	numControls = workspace.settings.optimalControlInfo[3]
	varShift = (numStates + numAlgebraic + numControls)*shiftSteps

	# check input
	@assert !(mode==:fullRTI) || workspace.settings.conservativismLevel >= 1
	@assert subsolverIterations >= -1
	@assert get_numVariables(newTailCosts) == varShift + numStates + numAlgebraic
	@assert get_numVariables(newTailConstraints) == varShift + numStates + numAlgebraic
	@assert mod(get_size(newTailVariables),numStates + numAlgebraic + numControls) == 0
	if newTailCosts.c != 0
		error("HBB mpc_shift!: The algorithm cannot handle objective function with constant parts, please reformulate the problem (possibly using an algebraic variable)")
	end

	# if there is no reference solution take as reference solution the best one in the tree
	if length(referenceSolution) == 0
        referenceSolutionNode = get_best_feasible_node(workspace,localOnly=localOnly)
        if referenceSolutionNode isa NullBBnode
            error("MPC shift: No suitable (partial) solution found for the current problem")
        else
            referenceSolution = referenceSolutionNode.primal
        end
    end

	# if there is no measured state take as measured state the state predicted by the reference solution (no model mismatch)
	if length(measuredState) == 0
		measuredState = referenceSolution[varShift+1:varShift+numStates]
	end

	@sync if !localOnly && !(workspace.sharedMemory isa NullSharedMemory)
         # call the local version of the function on the remote workers
         for p in workers()[1:workspace.settings.numProcesses-1]
             @async remotecall_fetch(OpenBB.eval,p,:(
			 	OpenBB.mpc_shift!(workspace,$shiftSteps,
								  $newTailCosts,$newTailConstraints,$newTailVariables,
								  referenceSolution=$referenceSolution,
								  measuredState=$measuredState,
								  mode=Symbol($(string(mode))),
								  approximateBellmanConstant=$approximateBellmanConstant,
								  subsolverIterations=$subsolverIterations,
								  suppressErrors=$suppressErrors,
                                  localOnly=true)
			 ))
         end
         # call the local version of the function on the main process
         mpc_shift!(workspace,shiftSteps,
		 		    newTailCosts,newTailConstraints,newTailVariables,
					referenceSolution=referenceSolution,
					measuredState=measuredState,
					mode=mode,
					approximateBellmanConstant=approximateBellmanConstant,
					subsolverIterations=subsolverIterations,
					suppressErrors=suppressErrors,
					localOnly=true)
	else

		# check for model mismatch (and adequateness of the settings)
		modelMismatch = false
		if maximum(abs.(referenceSolution[varShift+1:varShift+numStates] - measuredState)) > workspace.settings.primalTolerance
			@assert !(mode==:fullRTI) ||  workspace.settings.conservativismLevel == 2
			@assert !(mode==:fullRTI) || !(workspace.settings.withBoundsPropagation)
			modelMismatch = workspace.infeasiblesToRecover = true
		end

		# check for the necessity of recomputing the lower bounds
		if approximateBellmanConstant == Inf && mode != :constraintsOnly
			workspace.invalidLowerBounds = true
		end
		workspace.outdated = true

		# filter the tree depending on the chosen style
		if mode == :fullRTI
			# eliminate the nodes for which the selected discrete variables are not feasible
			dscShift = findfirst(ind->ind>varShift,workspace.problem.varSet.dscIndices)-1
			dscIndices_ = workspace.problem.varSet.dscIndices[1:dscShift]
			assignment = round.(referenceSolution[dscIndices_])
			function rti_selector(node::BBnode)::Bool
				for (k,ind) in enumerate(dscIndices_)
					if node.varLoBs[ind] > assignment[k] || node.varUpBs[ind] < assignment[k]
						return false
					end
				end
				return true
			end
			filter!(rti_selector,workspace.tree)
		elseif mode == :warmStart
			dscIndices = workspace.problem.varSet.dscIndices
			assignment = round.(referenceSolution[dscIndices])
			#ATTENTION: the following function has side effects (makes all nodes heuristic)
			function warmStart_selector!(node::BBnode)::Bool
				slack = 1.0
				for (k,ind) in enumerate(dscIndices)
					if node.varLoBs[ind] > assignment[k]
						slack -= node.varLoBs[ind] - assignment[k]
						if slack < 0.0 return false end
					elseif node.varUpBs[ind] < assignment[k]
						slack -= assignment[k] - node.varUpBs[ind]
						if slack < 0.0 return false end
					end
				end
				node.heuristic = true
				return true
			end
			filter!(warmStart_selector!,workspace.tree)
		elseif mode == :constraintsOnly
			# remove all nodes
			clear!(workspace)
		else
			throw(ArgumentError("BB mpc_shift!, shifting mode unknown"))
		end

		# make sure that the constraints are well sorted
		if !issorted_oc(workspace.problem.cnsSet,workspace.settings.optimalControlInfo)
			permutation = sortperm_oc(workspace.problem.cnsSet,workspace.settings.optimalControlInfo)
			permute_constraints!(workspace)
		end

		#### fix the initial part of the problem and impose measured state ####
		mpc_contraction!(workspace,shiftSteps,
						 measuredState,referenceSolution,
						 addTailVariables=is_empty(newTailVariables))


		#### Modify the tail of the problem #####
		_ABC = if mode!=:constraintsOnly approximateBellmanConstant else 0.0 end # it is useless to recompute bounds in :constraintsOnly mode
		mpc_extension!(workspace,shiftSteps,
					   referenceSolution,
					   newTailVariables=newTailVariables,
					   newTailConstraints=newTailConstraints,
					   newTailCosts=newTailCosts,
					   approximateBellmanConstant=_ABC)

		# update the workspace
		update!(workspace,updateAllNodes=true,subsolverIterations=subsolverIterations)

		# reinsert the root or update the nodes
		if mode == :fullRTI
			workspace.status = BBstatus(description="interrupted")
		elseif mode == :warmStart
			# reinsert the root
			push_node!(workspace,:active,BBroot(workspace))
			workspace.status = BBstatus(description="new")
		elseif mode == :constraintsOnly || mode == :warmStart
			# reset the updates register
			reset!(workspace.updatesRegister)
			# re-insert the root node into the queue
			if myid() == 1
				push_node!(workspace,:active,BBroot(workspace))
			end
			workspace.status = BBstatus(description="new")
		end

		# add a dummy node update to prevent wrong calculations of pseudoCosts
		push!(workspace.updatesRegister,do_nothing!,())

		if workspace.settings.verbose
			println("BB mpc_shift!, node count after shift."," Active: ",length(workspace.tree.active),",",
			                                                 " Suboptimals: ",length(workspace.tree.suboptimals),",",
															 " Infeasibles: ",length(workspace.tree.infeasibles),",",
															 " Solutions: ",length(workspace.tree.solutions),",",
															 " Blacklisted: ",length(workspace.tree.blacklisted),".")
		end
	end

	return
end


## Support functions

function mpc_contraction!(workspace::BBworkspace,shiftSteps::Int,measuredState::Vector{Float},referenceSolution::Vector{Float};addTailVariables::Bool=false)::Nothing

	# collect info on the problem
	numStates = workspace.settings.optimalControlInfo[1]
	varShift = sum(workspace.settings.optimalControlInfo)*shiftSteps
	dscIndices = copy(workspace.problem.varSet.dscIndices)

	# collect the first part of the objective
	if workspace.problem.objFun isa QuadraticObjective
		_Q = workspace.problem.objFun.Q[1:varShift,1:varShift]
		_L = workspace.problem.objFun.L[1:varShift]'
		objCorrection = x -> -.5*x'*_Q*x - _L*x
	else
		_L = workspace.problem.objFun.L[1:varShift]'
		objCorrection = x -> - _L*x
	end

	# fix the first part of the problem
	objectiveOffset = workspace.problem.objFun.c
	if !addTailVariables # remove the variable on the head of the problem (not preserving the definition for the new tail variables)
		fix_variables!(workspace.problem,1:varShift,referenceSolution[1:varShift],removeFixedVariables=true)
		tailVarsBounds = (Float[],Float[])
	else # shift backward the problem (preserving the definition for the new tail variables)
		shift_backward!(workspace.problem.cnsSet,varShift,referenceSolution[1:varShift])
		shift_backward!(workspace.problem.objFun,varShift,referenceSolution[1:varShift])
		tailVarsBounds = (workspace.problem.varSet.loBs[end-varShift+1:end],workspace.problem.varSet.upBs[end-varShift+1:end])
	end
	workspace.problem.objFun.c = objectiveOffset

	# remove the new constant constraints
	lastNZs = get_lastNZs(workspace.problem.cnsSet,2)
	toRemove = findall(iszero,lastNZs)
	remove_constraints!(workspace.problem.cnsSet,toRemove)

	# enforce the measured state
	workspace.problem.cnsSet.loBs[1:numStates] = measuredState
	workspace.problem.cnsSet.upBs[1:numStates] = measuredState

	# add the contraction function to the updates to perform on nodes
	push!(workspace.updatesRegister,contract!,(varShift,measuredState,referenceSolution,dscIndices,addTailVariables,tailVarsBounds,toRemove,objCorrection))

	return
end



function contract!(node::BBnode,varShift::Int,
				   measuredState::Vector{Float},
				   referenceSolution::Vector{Float},
				   dscIndices::Array{Int,1},
				   addTailVariables::Bool,
				   tailVarsBounds::Tuple{Vector{Float},Vector{Float}},
				   cnssToRemove::Vector{Int},
				   objCorrection::Function,
	      		   workspace::BBworkspace,
				   newUpdatedVars::BitArray{1},
				   newUpdatedCnss::BitArray{1})::Bool

	# collect info
	numStates = length(measuredState)

	#### Constraints ####
	deleteat!(node.cnsLoBs,cnssToRemove)
	deleteat!(node.cnsUpBs,cnssToRemove)
	deleteat!(node.cnsDual,cnssToRemove)
	node.cnsLoBs[1:numStates] = measuredState
	node.cnsUpBs[1:numStates] = measuredState


	#### Variables ####
	# reduce fractionality (due to variable elimination)
	dscShift = findfirst(ind->ind>varShift,dscIndices)-1
	deltaAF = @. abs(node.primal[dscIndices[1:dscShift]] - round(node.primal[dscIndices[1:dscShift]]))
	@. deltaAF =  deltaAF*(deltaAF>workspace.settings.primalTolerance)
	node.avgFractionality -=  2.0*sum(deltaAF)/length(dscIndices)

	# reduce pseudo cost (due to variable elimination)
    for (k,i) in enumerate(dscIndices[1:dscShift])
        node.pseudoCost -= min(workspace.problem.varSet.pseudoCosts[1][k,1]*(node.primal[i]-floor(node.primal[i]+workspace.settings.primalTolerance)),
                               workspace.problem.varSet.pseudoCosts[1][k,2]*(ceil(node.primal[i]-workspace.settings.primalTolerance)-node.primal[i]))
    end

	# perform the shift
	if !addTailVariables

		# remove head from the primal solution
		deleteat!(node.primal,1:varShift)

		# remove head from the dual solution
		deleteat!(node.bndDual,1:varShift)

		# remove head from variable bounds
		deleteat!(node.varLoBs,1:varShift)
		deleteat!(node.varUpBs,1:varShift)

		# remove the variables from the cuts
		fix_variables!(node.cutSet,1:varShift,referenceSolution[1:varShift],removeFixedVariables=true)

	else

		# shift backward the primal solution
		oldPrimal = node.primal[1:varShift+numStates]
		node.primal[1:end-varShift] = node.primal[varShift+1:end]
		node.primal[1:numStates] = measuredState

		# shift backward the variable bounds
		node.varLoBs[1:end-varShift] = node.varLoBs[varShift+1:end];
		node.varLoBs[end-varShift+1:end] = tailVarsBounds[1]
		node.varUpBs[1:end-varShift] = node.varUpBs[varShift+1:end]
		node.varUpBs[end-varShift+1:end] = tailVarsBounds[2]

		# shift backward the bound-duals
		node.bndDual[1:end-varShift] = node.bndDual[varShift+1:end]
		node.bndDual[end-varShift+1:end] .= 0.0

		# shift backward the cuts
		if get_size(node.cutSet) > 0
			shift_backward!(node.cutSet,varShift,referenceSolution)
		end

		# correct the objective lower bound
		node.objLoB += objCorrection(oldPrimal[1:varShift]) + node.cnsDual[1:numStates]'*(oldPrimal[varShift+1:end]-measuredState)
		node.objUpB = Inf

		# assume the new tail variables to be maximally fractional
		dscTailLimit = findlast(ind->ind<length(node.primal)-varShift,dscIndices)+1
		node.avgFractionality +=  (dscTailLimit-length(dscIndices))/length(dscIndices)
		for k in dscTailLimit:length(dscIndices)
			node.pseudoCost += 0.5*min(workspace.problem.varSet.pseudoCosts[1][k,1],workspace.problem.varSet.pseudoCosts[1][k,2])
		end
	end

	# remove the null cuts
	if get_size(node.cutSet) > 0
		toRemove = findall(iszero,get_lastNZs(node.cutSet,2))
		if !isempty(toRemove)
			remove_constraints!(node.cutSet,toRemove)
			deleteat!(node.cutDual,toRemove)
		end
	end

	#### Preprocessing #### (although it shouldn't happen at this stage)
	# reove the removed constraints
	deleteat!(newUpdatedCnss,cnssToRemove)
	newUpdatedCnss[1:numStates] .= true

	# remove/shift the variables and constraints to preprocess
    if !addTailVariables
		deleteat!(newUpdatedVars,1:varShift)
	else
		newUpdatedVars[numStates+1:end-varShift] = newUpdatedVars[varShift+numStates+1:end]
		newUpdatedVars[end-varShift+1:end] .= true
	end

	return true
end



function mpc_extension!(workspace::BBworkspace,extensionSteps::Int,
						referenceSolution::Array{Float64,1};
						newTailVariables::VariableSet=EmptyVariableSet(),
						newTailConstraints::ConstraintSet=NullConstraintSet(),
						newTailCosts::ObjectiveFunction=NullConstraintSet(),
						approximateBellmanConstant)::Nothing

	# collect info
	numVars = get_size(workspace.problem.varSet) + get_size(newTailVariables)
	(numStates,numAlgebraic,numControls) = workspace.settings.optimalControlInfo
	numVarsPerStep = (numStates + numAlgebraic + numControls)

	oldNumOfSteps = div((numControls-1)+maximum(get_lastNZs(workspace.problem.cnsSet,2)),numVarsPerStep)
	newNumOfSteps = oldNumOfSteps+extensionSteps

	rangeOldTerminalVars = oldNumOfSteps*numVarsPerStep+1:oldNumOfSteps*numVarsPerStep+numStates+numAlgebraic 	# current position of the terminal objective terms
	rangeNewTerminalVars = newNumOfSteps*numVarsPerStep+1:newNumOfSteps*numVarsPerStep+numStates+numAlgebraic 	# future position of the terminal terms

	oldTimeSteps = get_timeSteps_oc(workspace.problem.cnsSet,(numStates,numAlgebraic,numControls))
	rangeOldTerminalCnss = findlast((x)->(x<=oldNumOfSteps),oldTimeSteps)+1:get_size(workspace.problem.cnsSet)

	locationNewVars = oldNumOfSteps*numVarsPerStep+1:newNumOfSteps*numVarsPerStep+numStates+numAlgebraic      # future position of the new steps + terminal terms


	# check size consistency
	@assert newNumOfSteps*numVarsPerStep + numStates + numAlgebraic == numVars
	@assert get_numVariables(newTailConstraints) == extensionSteps*numVarsPerStep + numStates + numAlgebraic
	@assert get_numVariables(newTailCosts) == extensionSteps*numVarsPerStep + numStates + numAlgebraic

   # check type consistency
	if workspace.problem.objFun isa LinearObjective
		@assert is_linear(newTailCosts)
	elseif workspace.problem.objFun isa QuadraticObjective
		@assert is_quadratic(newTailCosts)
	end

	# check for quadratic parts in the terminal objective (not allowed)
	if newTailCosts isa QuadraticObjective && !iszero(newTailCosts.Q[end-numStates-numAlgebraic+1:end,end-numStates-numAlgebraic+1:end])
		throw(ArgumentError("There cannot be a quadratic part in the terminal objective"))
	end


	################## Variables ##################
	# extend the problem if needed
	numNewVariables = get_size(newTailVariables)
	if numNewVariables > 0
		newTailVarsLoBs = newTailVariables.loBs
		newTailVarsUpBs = newTailVariables.upBs
		dscCounts = (get_numDiscrete(workspace.problem.varSet),get_numDiscrete(newTailVariables))
		append!(workspace.problem.varSet,newTailVariables)
		append_variables!(workspace.problem.cnsSet,numNewVariables)
		append_variables!(workspace.problem.objFun,numNewVariables)
	else
		newTailVarsLoBs = Float[]
		newTailVarsUpBs = Float[]
		dscCounts = (get_numDiscrete(workspace.problem.varSet),0)
	end

	################## Constraints ##################

	# store the terminal constraints
	oldTerminalCnss = workspace.problem.cnsSet[rangeOldTerminalCnss]

	# remove the old terminal constraints
	remove_constraints!(workspace.problem.cnsSet,rangeOldTerminalCnss)

    # append the new constraints to the constraint set
	newTailConstraints_ = deepcopy(newTailConstraints)
	insert_variables!(newTailConstraints_,oldNumOfSteps*numVarsPerStep,1)

	# update the constraint set
    append!(workspace.problem.cnsSet,newTailConstraints_)


	################## Objective ##################

	# collect the objective terms to add
	newTailCosts_ = deepcopy(newTailCosts)
	insert_variables!(newTailCosts_,oldNumOfSteps*numVarsPerStep,1)

	# update the objective
	workspace.problem.objFun.L[locationNewVars] = newTailCosts_.L[locationNewVars]
	if workspace.problem.objFun isa QuadraticObjective
		workspace.problem.objFun.Q[locationNewVars,locationNewVars] = newTailCosts_.Q[locationNewVars,locationNewVars]
		workspace.problem.objFun.pInvQ = pinv(workspace.problem.objFun.Q)
	end


	############## Nodes Update (with dual correction) #####################
	if approximateBellmanConstant == Inf
	    # create an optimization problem to compute a correction term for the nodes duals
		numCnss = get_size(workspace.problem.cnsSet)
		correctionProblem = Problem(varSet=deepcopy(workspace.problem.varSet),
									cnsSet=deepcopy(workspace.problem.cnsSet),
									objFun=newTailCosts_)

		# build a temporary node as results space this optimization
		node_ = BBroot(workspace)

		# initialize the matrices that will be used to store the dual correction terms
		bndCorrectionMatrix = Array{Float64,2}(undef,numVars,max(1,length(rangeOldTerminalCnss))) # bounds
		cnsCorrectionMatrix = Array{Float64,2}(undef,numCnss,max(1,length(rangeOldTerminalCnss))) # constraints

		if isempty(rangeOldTerminalCnss)
			# solve correction problem as it is and store the results
			out = solve!(node_,setup(correctionProblem,workspace.subsolverWS.settings))
			bndCorrectionMatrix[:,1] = node_.bndDual
			cnsCorrectionMatrix[:,1] = node_.cnsDual

		else

			# build and solve the correction problems
			for k in 1:length(rangeOldTerminalCnss)

				# add one of the old terminal constraints as a cost
				correctionProblem.objFun.L[rangeOldTerminalVars[1:numStates]] += -oldTerminalCnss.A[k,rangeOldTerminalVars[1:numStates]]

				# solve correction problem and store the results
				out = solve!(node_,setup(correctionProblem,workspace.subsolverWS.settings))
				if out[1] == 0
					bndCorrectionMatrix[:,k] = node_.bndDual
					cnsCorrectionMatrix[:,k] = node_.cnsDual
				else
					@warn "Issue in the dual recomputation encountered"
					bndCorrectionMatrix[:,k] .= 0.0
					cnsCorrectionMatrix[:,k] .= 0.0
				end
			end
		end
	else
		bndCorrectionMatrix = Array{Float64,2}(undef,0,0)
		cnsCorrectionMatrix = Array{Float64,2}(undef,0,0) # constraints
	end

	# add the extension function to the updates to perform on nodes
	push!(workspace.updatesRegister,extend!,(rangeOldTerminalCnss,
											 bndCorrectionMatrix,cnsCorrectionMatrix,
											 newTailVarsLoBs,newTailVarsUpBs,dscCounts,
											 newTailConstraints.loBs,newTailConstraints.upBs,
											 approximateBellmanConstant))


	# sort again the constraints
	sort_constraints_oc!(workspace,localOnly=true)

	return

end


function extend!(node::BBnode,rangeOldTerminalCnss::UnitRange{Int64},
				 bndCorrectionMatrix::Array{Float64,2},cnsCorrectionMatrix::Array{Float64,2},
				 newVarsLoBs::Array{Float64,1},newVarsUpBs::Array{Float64,1},dscCounts::Tuple{Int,Int},
				 newCnssLoBs::Array{Float64,1},newCnssUpBs::Array{Float64,1},approximateBellmanConstant::Float,
				 workspace::BBworkspace,newUpdatedVars::BitArray{1},newUpdatedCnss::BitArray{1})::Bool

	########## Variables #######
	if !isempty(newVarsLoBs)
		append!(node.varLoBs,newVarsLoBs)
		append!(node.varUpBs,newVarsUpBs)
		append!(node.primal,zeros(length(newVarsLoBs)))
		append!(node.bndDual,zeros(length(newVarsLoBs)))
		append!(newUpdatedVars,ones(Bool,length(newVarsLoBs)))
	end

	# adapt fractionality
 	node.avgFractionality = (node.avgFractionality*dscCounts[1] + dscCounts[2])/sum(dscCounts)


	######## Constraints #######
	# remove the duals relative to the terminal constraints
	terminalCnsDuals = splice!(node.cnsDual,rangeOldTerminalCnss)
	deleteat!(node.cnsLoBs,rangeOldTerminalCnss)
	deleteat!(node.cnsUpBs,rangeOldTerminalCnss)

	# expand the constraint info
	append!(node.cnsDual,zeros(length(newCnssLoBs)))
	append!(node.cnsLoBs,copy(newCnssLoBs))
	append!(node.cnsUpBs,copy(newCnssUpBs))

	# compute the correction for the duals
	if approximateBellmanConstant == Inf
		terminalCnsDuals = terminalCnsDuals./sum(terminalCnsDuals)
		if isempty(rangeOldTerminalCnss)
			node.bndDual .+= dropdims(bndCorrectionMatrix,dims=2)
			node.cnsDual .+= dropdims(cnsCorrectionMatrix,dims=2)
		else
			node.bndDual .+= bndCorrectionMatrix*terminalCnsDuals
			node.cnsDual .+= cnsCorrectionMatrix*terminalCnsDuals
		end
	else
		node.objLoB -= approximateBellmanConstant
	end

	return true
end
