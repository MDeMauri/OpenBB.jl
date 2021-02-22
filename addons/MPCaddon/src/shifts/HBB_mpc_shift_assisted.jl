# @Author: Massimo De Mauri <massimo>
# @Date:   2021-01-15T21:26:18+01:00
# @Email:  massimo.demauri@protonmail.com
# @Filename: HBB_mpc_shift_assisted.jl
# @Last modified by:   massimo
# @Last modified time: 2021-02-22T17:08:55+01:00
# @License: LGPL-3.0
# @Copyright: {{copyright}}


# this version of the mpc shift expects as an already shifted MINLP and performs
# an actual shift only on the linear relaxation already built. Depending on
# how the shifted problem is obtained, this function can be more efficient than
# the actual non-linear shifting procedure. However, it is a dangerous function
# it is possible to mess the definition of the shifted problem.

function _mpc_shift!(workspace::HBBworkspace,shiftSteps::Int,shiftedMINLP::Problem;
					referenceSolution::Array{Float64,1}=Float64[],measuredState::Array{Float64,1}=Float64[],
					mode::Symbol=:fullRTI,approximateBellmanConstant::Float=Inf,subsolverIterations::Int=0,
	                suppressErrors::Bool=false)::Nothing

	numVars = get_numVariables(workspace)
	(numStates,numAlgebraic,numControls) = workspace.settings.optimalControlInfo
	numVarsPerStep = numStates+numAlgebraic+numControls
	varShift = shiftSteps*numVarsPerStep
	lastNZs = get_lastNZs(workspace.problem.cnsSet,2)
	numSteps = div((numControls-1)+maximum(lastNZs),numVarsPerStep)
	timeSteps = get_timeSteps_oc(shiftedMINLP.cnsSet,(numStates,numAlgebraic,numControls))
	tailCnssInds = findlast(x->x<=numSteps-shiftSteps,timeSteps)+1:length(timeSteps)
	tailVarsInds = (numSteps-shiftSteps)*numVarsPerStep+1:numVars

	# collect tail variables
	if workspace.problem.varSet != shiftedMINLP.varSet
		tailVariables = deepcopy(shiftedMINLP.varSet)
		remove_variables!(tailVariables,1:numVars-varShift)
	else
		tailVariables = EmptyVariableSet()
	end

	 # collect the tail costs
	tailCosts = deepcopy(shiftedMINLP.objFun)
	remove_variables!(tailCosts,1:numVars-varShift-numStates-numAlgebraic)

	# collect the tail constraints
	tailConstraints = shiftedMINLP.cnsSet[tailCnssInds]
	remove_variables!(tailConstraints,1:numVars-varShift-numStates-numAlgebraic)

	# call normal shift
	mpc_shift!(workspace,shiftSteps,tailCosts,tailConstraints,tailVariables;
			   referenceSolution=referenceSolution,measuredState=measuredState,
			   mode=mode,approximateBellmanConstant=approximateBellmanConstant,
			   subsolverIterations=subsolverIterations,suppressErrors=suppressErrors)::Nothing
	return
end



function mpc_shift!(workspace::HBBworkspace,shiftSteps::Int,shiftedMINLP::Problem;
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

	# if there is no reference solution assume as
	# reference solution the best one in the tree
	if length(referenceSolution) == 0
		if !isempty(workspace.feasibleNodes)
			referenceSolution = workspace.feasibleNodes[1].primal
		else
			#TODO check if we can borrow some partial solution from BB
			error("HBB, mpc_shift!: No suitable (partial) solution found for the current problem")
		end
	end

	# if there is no new initial state assume as new initial state
	# the state predicted by the reference solution
	if isempty(measuredState)
		measuredState = referenceSolution[varShift+1:varShift+numStates]
	end

	# check for model mismatch (and adequateness of the settings)
	inexactModel = maximum(abs.(referenceSolution[varShift+1:varShift+numStates] - measuredState)) > workspace.settings.primalTolerance
	if inexactModel
		@assert workspace.settings.conservativismLevel == 2
	end

	# make sure that the constraints are well sorted
	permutation = sortperm_oc(shiftedMINLP.cnsSet,workspace.settings.optimalControlInfo)
	if !issorted(permutation)
		permute!(shiftedMINLP.cnsSet,localOnly=true)
	end


	# collect info on the tail costs and constraints
	lastNZs = get_lastNZs(workspace.problem.cnsSet,2)
	numSteps = div((numControls-1)+maximum(lastNZs),numVarsPerStep)
	timeSteps = get_timeSteps_oc(shiftedMINLP.cnsSet,(numStates,numAlgebraic,numControls))
	tailCnssInds = findlast(x->x<=numSteps-shiftSteps,timeSteps)+1:length(timeSteps)
	tailVarsInds = (numSteps-shiftSteps)*numVarsPerStep+1:numVars

 	# collect the tail costs
	tailCosts = deepcopy(shiftedMINLP.objFun)
	remove_variables!(tailCosts,1:(numSteps-shiftSteps)*numVarsPerStep)

	# collect the tail constraints
	tailConstraints = shiftedMINLP.cnsSet[tailCnssInds]
	remove_variables!(tailConstraints,1:(numSteps-shiftSteps)*numVarsPerStep)

	# collect info for linear relaxation of tail constraints
	linearInds = findall(iszero,tailConstraints.hesSparsity)
	nonLinearInds = findall(!iszero,tailConstraints.hesSparsity)

	# extract linear tail constraints
	tailRelaxations = linearRelaxation(tailConstraints,zeros(length(tailVarsInds)))[linearInds]
	relaxationsIDs = tailCnssInds[linearInds].+numCnss

	# generate linear relaxations for the non-linear tail constraints
	nonLinearTailCnss = tailConstraints[nonLinearInds]
	terminalNonLinearInds = findall(ts->ts>shiftSteps,get_timeSteps_oc(nonLinearTailCnss,workspace.settings.optimalControlInfo))
	tolerance = workspace.settings.primalTolerance
	for node in workspace.feasibleNodes
		cnssVal = evaluate(nonLinearTailCnss,node.primal[tailVarsInds])
		cnsActivation = @. !(nonLinearTailCnss.loBs+tolerance <= cnssVal <= nonLinearTailCnss.upBs-tolerance)
		active = union(findall(!iszero,cnsActivation),terminalNonLinearInds)
		append!(tailRelaxations,linearRelaxation(nonLinearTailCnss,node.primal[tailVarsInds])[active])
		append!(relaxationsIDs,tailCnssInds[nonLinearInds[active]].+numCnss)
	end

	# sort the new linearizations
	permutation = sortperm_oc(tailRelaxations,workspace.settings.optimalControlInfo)
	if !issorted(permutation)
		permute!(tailRelaxations,permutation)
		permute!(relaxationsIDs,permutation)
	end

	# update the workspace
	empty!(workspace.feasibleNodes)
	workspace.status = HBBstatus(description="interrupted")

	# propagate the shift to the mixed-integer subsolver
	mpc_shift!(workspace.mipStepWS,shiftSteps,tailCosts,tailRelaxations,relaxationsIDs,
			   referenceSolution=referenceSolution,measuredState=measuredState,mode=mode,
			   approximateBellmanConstant=approximateBellmanConstant,
			   subsolverIterations=subsolverIterations,
			   suppressErrors=suppressErrors)


	# constraint removal (head and tail)
	toRemove = vcat(findall(x->x<=numVarsPerStep*shiftSteps,lastNZs),findall(x->x>numSteps,timeSteps))
	# here in theory: remove_constraints!(workspace,toRemove); but instead:
	oldNumRelaxations = get_numConstraints(workspace.mipStepWS.mipSolverWS) #TOREMOVE: debug
	remove_constraints!(workspace.mipStepWS,toRemove,suppressErrors=true)
	@assert get_numConstraints(workspace.mipStepWS.mipSolverWS) == oldNumRelaxations #TOREMOVE: debug

	# here in theory we have to shift backward the problem and then add the new tail stuff to the problem
	# but instead: substitute the new problem to the old
	@assert get_numVariables(shiftedMINLP) == get_numVariables(workspace.problem)
	workspace.problem.varSet = deepcopy(shiftedMINLP.varSet)
	workspace.problem.cnsSet = deepcopy(shiftedMINLP.cnsSet)
	workspace.problem.objFun = deepcopy(shiftedMINLP.objFun)

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
