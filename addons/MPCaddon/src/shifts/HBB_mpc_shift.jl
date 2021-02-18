# @Author: Massimo De Mauri <massimo>
# @Date:   2021-01-16T16:18:21+01:00
# @Email:  massimo.demauri@protonmail.com
# @Filename: HBB_mpc_shift.jl
# @Last modified by:   massimo
# @Last modified time: 2021-02-16T19:29:32+01:00
# @License: LGPL-3.0
# @Copyright: {{copyright}}


function mpc_shift!(workspace::HBBworkspace,shiftSteps::Int,newTailCosts::ObjectiveFunction,
					newTailConstraints::ConstraintSet,newTailVariables::VariableSet=EmptyVariableSet(),shiftedMINLP=nothing;
					referenceSolution::Array{Float64,1}=Float64[],measuredState::Array{Float64,1}=Float64[],
					mode::Symbol=:fullRTI,approximateBellmanConstant::Float=Inf,subsolverIterations::Int=0,
					suppressErrors::Bool=false)::Nothing

	# check workspace
	if all(workspace.settings.optimalControlInfo .== -1)
		error("mpc_shift!: The HBBworkspace has no optimal control info.")
	else
		@assert all(workspace.settings.optimalControlInfo .>= 0)
	end

	# check the problem
	if workspace.problem.objFun.c != 0
		error("HBB mpc_shift!: The algorithm cannot handle objective function with constant parts, please reformulate the problem (possibly using an algebraic variable)")
	end

	# collect info
	(numStates,numAlgebraic,numControls) = workspace.settings.optimalControlInfo
	numVarsPerStep = numStates + numAlgebraic + numControls
	varShift = numVarsPerStep*shiftSteps
	numVars = get_numVariables(workspace.problem)
	numCnss = get_numConstraints(workspace.problem)

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
		if !isempty(workspace.feasibleNodes)
			referenceSolution = workspace.feasibleNodes[1].primal
		else
			#TODO check if we can borrow some partial solution from BB
			error("HBB, mpc_shift!: No suitable (partial) solution found for the current problem")
		end
	end

	# if there is no new initial state take as new initial state the state predicted by the reference solution
	if isempty(measuredState)
		measuredState = referenceSolution[varShift+1:varShift+numStates]
	end

	# check for model mismatch (and adequateness of the settings)
	deltaState = measuredState - referenceSolution[varShift+1:varShift+numStates]
	@assert maximum(abs.(deltaState)) < workspace.settings.primalTolerance || workspace.settings.conservativismLevel == 2

	# make sure that the constraints are well sorted
	_permutation = sortperm_oc(workspace.problem.cnsSet,workspace.settings.optimalControlInfo)
	if !issorted(_permutation) permute_constraints!(workspace,_permutation) end
	_permutation = sortperm_oc(newTailConstraints,workspace.settings.optimalControlInfo)
	if !issorted(_permutation) permute!(newTailConstraints,_permutation) end

	# collect info for linear relaxation of tail constraints
	linearInds = findall(iszero,newTailConstraints.hesSparsity); linearIDs = linearInds .+ numCnss
	nonLinearInds = findall(!iszero,newTailConstraints.hesSparsity); nonLinearIDs = linearInds .+ numCnss

	# extract linear tail constraints
	newTailRelaxations = linearRelaxation(newTailConstraints,zeros(varShift+numStates+numAlgebraic))[linearInds]
	relaxationsIDs = copy(linearIDs)

	# generate linear relaxations for the non-linear tail constraints
	nonLinearTailCnss = newTailConstraints[nonLinearInds]
	terminalNonLinearInds = findall(ts->ts>shiftSteps,get_timeSteps_oc(nonLinearTailCnss,workspace.settings.optimalControlInfo))
	tolerance = workspace.settings.primalTolerance
	for node in workspace.feasibleNodes
		cnssVal = evaluate(nonLinearTailCnss,node.primal[numVars-varShift-numStates-numAlgebraic+1:numVars])
		_activeMask = @. !(nonLinearTailCnss.loBs+tolerance <= cnssVal <= nonLinearTailCnss.upBs-tolerance)
		_cnsSelection = union(findall(!iszero,_activeMask),terminalNonLinearInds)
		append!(newTailRelaxations,linearRelaxation(nonLinearTailCnss,node.primal[numVars-varShift-numStates-numAlgebraic+1:numVars])[_cnsSelection])
		append!(relaxationsIDs,nonLinearIDs[_cnsSelection])
	end

	# sort the new linearizations
	_permutation = sortperm_oc(newTailRelaxations,workspace.settings.optimalControlInfo)
	if !issorted(_permutation)
		permute!(newTailRelaxations,_permutation)
		permute!(relaxationsIDs,_permutation)
	end

	# propagate the shift to the mixed-integer subsolver
	mpc_shift!(workspace.mipStepWS,shiftSteps,newTailCosts,newTailRelaxations,relaxationsIDs,
			   newTailVariables,referenceSolution=referenceSolution,
			   measuredState=measuredState,mode=mode,
			   approximateBellmanConstant=approximateBellmanConstant,
			   subsolverIterations=subsolverIterations,
			   suppressErrors=suppressErrors)



	# collect info on the problem costs and constraints
	lastNZs = get_lastNZs(workspace.problem.cnsSet,2)
	numSteps = div((numControls-1)+maximum(lastNZs),numVarsPerStep)
	timeSteps = get_timeSteps_oc(workspace.problem.cnsSet,(numStates,numAlgebraic,numControls))

	# shift variables
	if get_size(newTailVariables) > 0
		remove_variables!(workspace.problem.varSet,1:varShift)
		append!(workspace.problem.varSet,newTailVariables)
	end

	# shift constraints
	_oldHeadTailCnssInds = vcat(findall(x->x<=numVarsPerStep*shiftSteps,lastNZs),findall(x->x>numSteps,timeSteps))
	# debug
	# remove_constraints!(workspace.mipStepWS,_oldHeadTailCnssInds,suppressErrors=true)
	# workspace.problem.cnsSet = deepcopy(shiftedMINLP.cnsSet)
	# debug
	fix_variables!(workspace.problem.cnsSet,1:varShift,referenceSolution[1:varShift],removeFixedVariables=true) # remove head variables
	append_variables!(workspace.problem.cnsSet,varShift) # append new variables at the end
	_newTailCnss = deepcopy(newTailConstraints); insert_variables!(_newTailCnss,numVars-varShift-numStates-numAlgebraic,1) # adapt the size of the new tail constraints
	append!(workspace.problem.cnsSet,_newTailCnss) # append the new tail constraints
	remove_constraints!(workspace,_oldHeadTailCnssInds,suppressErrors=true) # old remove tail/head constraints

	# impose measured state
	_shiftedSolution = copy(referenceSolution); _shiftedSolution[1:end-varShift] = _shiftedSolution[varShift+1:end]; _shiftedSolution[1:numStates] = measuredState
	workspace.problem.cnsSet.loBs[1:numStates] =
	workspace.problem.cnsSet.upBs[1:numStates] = evaluate(workspace.problem.cnsSet,_shiftedSolution)[1:numStates]

	@info linearRelaxation(workspace.problem.cnsSet,_shiftedSolution)[1:numStates]

	# shift objective
	fix_variables!(workspace.problem.objFun,numVars-numStates-numAlgebraic+1:numVars,zeros(numStates+numAlgebraic),removeFixedVariables=true) # remove tail variables
	fix_variables!(workspace.problem.objFun,1:varShift,referenceSolution[1:varShift],removeFixedVariables=true); workspace.problem.objFun.c = 0 # remove head variables
	append_variables!(workspace.problem.objFun,varShift+numStates+numAlgebraic) # append new variables at the end
 	_newTailCosts = deepcopy(newTailCosts); insert_variables!(_newTailCosts,numVars-varShift-numStates-numAlgebraic,1) # adapt the size of the new tail costs
  	add!(workspace.problem.objFun,_newTailCosts) # add the new tail costs

	# update the workspace
	empty!(workspace.feasibleNodes)
	workspace.linearCnssIndices = findall(iszero,workspace.problem.cnsSet.hesSparsity)
	workspace.nonLinearCnssIndices = collect(1:get_size(workspace.problem.cnsSet))
	deleteat!(workspace.nonLinearCnssIndices,workspace.linearCnssIndices)
	workspace.status = HBBstatus(description="interrupted")
	workspace.blackList = BBblackList(get_numDiscrete(workspace.problem))
	workspace.outdated = false

	# mark the remaining workspace components as outdated
	make_outdated!(workspace.nlpStepWS)
	make_outdated!(workspace.heuristicsWS)

	return
end


## Support Functions
# MPC shift function for the mip step
function mpc_shift!(workspace::HBBmipStepWorkspace,shiftSteps::Int,newTailCosts::ObjectiveFunction,newTailConstraints::ConstraintSet,newTailConstraintsIDs::Vector{Int},
					newTailVariables::VariableSet=EmptyVariableSet();referenceSolution::Array{Float64,1}=Float64[],measuredState::Array{Float64,1}=Float64[],
					mode::Symbol=:fullRTI,approximateBellmanConstant::Float=Inf,subsolverIterations::Int=0,suppressErrors::Bool=false)::Nothing

	# collect info
	(numStates,numAlgebraic,numControls) = workspace.mipSolverWS.settings.optimalControlInfo
	numVarsPerStep = sum(workspace.mipSolverWS.settings.optimalControlInfo)

	# make sure that the relaxations are well sorted
	_permutation = sortperm_oc(workspace.mipSolverWS.problem.cnsSet,workspace.mipSolverWS.settings.optimalControlInfo)
	if !issorted(_permutation)
		permute_constraints!(workspace.mipSolverWS,_permutation)
		permute!(workspace.relaxationsIDs,_permutation)
	end
	_permutation = sortperm_oc(newTailConstraints,workspace.mipSolverWS.settings.optimalControlInfo)
	if !issorted(_permutation)
		permute_constraints!(newTailConstraints,_permutation)
		permute!(newTailConstraintsIDs,_permutation)
	end

	# remove the IDs of the relaxations that will become constant and of the terminal constraints
	lastNZs = get_lastNZs(workspace.mipSolverWS.problem.cnsSet,2)
	numSteps = div((numControls-1)+maximum(lastNZs),numVarsPerStep)
	timeSteps = get_timeSteps_oc(workspace.mipSolverWS.problem.cnsSet,(numStates,numAlgebraic,numControls))
	toRemove = vcat(findall(x->x<=numVarsPerStep*shiftSteps,lastNZs),findall(x->x>numSteps,timeSteps))
	deleteat!(workspace.relaxationsIDs,toRemove)

	# the new relaxations got on top of all the rest
	append!(workspace.relaxationsIDs,newTailConstraintsIDs)

	# clear the parameters set by HBB on BB
	workspace.mipSolverWS.settings.objectiveCutoff = Inf
	clear_blackList!(workspace.mipSolverWS)

	# actually perform the MPC shift
	mpc_shift!(workspace.mipSolverWS,shiftSteps,newTailCosts,newTailConstraints,newTailVariables,
			   referenceSolution=referenceSolution,measuredState=measuredState,mode=mode,
			   approximateBellmanConstant=approximateBellmanConstant,subsolverIterations=subsolverIterations,
			   suppressErrors=suppressErrors)

	return
end
