# @Author: Massimo De Mauri <massimo>
# @Date:   2020-11-25T13:23:16+01:00
# @Email:  massimo.demauri@gmail.com
# @Filename: solve.jl
# @Last modified by:   massimo
# @Last modified time: 2021-01-28T16:31:10+01:00
# @License: LGPL-3.0
# @Copyright: {{copyright}}


function solve!(workspace::HBBworkspace,initialPoint::Vector{Float}=Float[])::Nothing

	# set the workspace in running state and start collecting time
	initialTime = time()
	workspace.status.description = "running"
	if workspace.settings.verbose println("HBB: Started...") end

	# synchronize the libraries with the remote workers if necessary
	numProcesses = max(workspace.settings.nlpProcesses,workspace.settings.mipProcesses)
	if numProcesses > 1 && libsManager.modified
		libsManager.modified=false
		@sync for k in 2:numProcesses
			@async remotecall_fetch(OpenBB.eval,k,:(OpenBB.remove_all_libs();OpenBB.load_libs($(list_libs()))))
		end
	end


	# get info on the problem
	numVars = get_size(workspace.problem.varSet)
	numCnss	= get_size(workspace.problem.cnsSet)
	dscIndices = workspace.problem.varSet.dscIndices
	(cnsLoBs,cnsUpBs) = get_bounds(workspace.problem.cnsSet)
	cnsJacobianType = get_jacobianType(workspace.problem.cnsSet)


	###################  update the workspace if necessary ##################
	if workspace.outdated
		guessList = workspace.feasibleNodes
		update!(workspace)
	else
		guessList = BBnode[]
	end

	################## generate an initial point if needed ######################
	if isempty(guessList)
		if isempty(workspace.feasibleNodes) && isempty(initialPoint)

			# generate dummy initial node
			guessList = [BBnode(copy(workspace.problem.varSet.loBs),copy(workspace.problem.varSet.upBs),
							    copy(workspace.problem.cnsSet.loBs),copy(workspace.problem.cnsSet.upBs),
								zeros(numVars),zeros(numVars),zeros(numCnss))]

			# solve nlp problem to get a better initial node
			nlpWS_ = setup(workspace.problem,workspace.settings.nlpSettings)
			solve!(guessList[1],nlpWS_)

		elseif !isempty(initialPoint)
			@assert length(initialPoint) == numVars
			guessList = [BBnode(copy(workspace.problem.varSet.loBs),copy(workspace.problem.varSet.upBs),
								copy(workspace.problem.cnsSet.loBs),copy(workspace.problem.cnsSet.upBs),
								initialPoint,zeros(numVars),zeros(numCnss))]
		else
			guessList = [deepcopy(workspace.feasibleNodes[k]) for k in 1:min(length(workspace.feasibleNodes),workspace.settings.nlpProcesses)]
		end
	end

	############################# main iteration ####################
	numIterations = 0
	cnsJacobian = zeros(numCnss,numVars)
	cnsJacSparsity = get_jacobianSparsity(workspace.problem.cnsSet)
	repeatedDetector = node->lookup(workspace.blackList,node.primal[dscIndices])
	while true

		# update counters
		numIterations += 1; elapsedTime = time()-initialTime

		# print status
		if workspace.settings.verbose
			println("HBB iter "*string(numIterations)*":"*
			" Obj.Bounds=("*string(workspace.status.objLoB)*","*string(workspace.status.objUpB)*")"*
			" Abs.Gap="*string(workspace.status.absoluteGap)*" Rel.Gap="*string(round(workspace.status.relativeGap,digits=2)))
	   	end

		######################### NLP ########################
		timeCheckpoint = time()
		if workspace.settings.verbose println("- NLP step:") end

		if !isempty(guessList)

			# add the current guesses to the black-list to check if they are found more than once
			insert!(workspace.blackList,getindex.(getfield.(guessList,:primal),[dscIndices]))

			# perform nlp step
			for guess in guessList
				# construct a guess node of appropriate dimensions and containing a guess for the constraint duals
				jacc_ = evaluate_jacobian(workspace.problem.cnsSet,guess.primal)
				if jacc_ isa SpMatrix{Float} # pinv works only on sparse matrices
					@. cnsJacobian[cnsJacSparsity[1],cnsJacSparsity[2]] = jacc_[cnsJacSparsity[1],cnsJacSparsity[2]]
				else
					cnsJacobian = jacc_
				end
				guess.cnsDual = pinv(cnsJacobian')*(-evaluate_gradient(workspace.problem.objFun,guess.primal)-guess.bndDual)
				(guess.cnsLoBs,guess.cnsUpBs) = get_bounds(workspace.problem.cnsSet)

				# # add a cut to the node in order to enforce a maximum level for the objective
				# if guess.maxNumOfCuts > get_size(guess.cutSet)
				# 	objApproximation = firstOrderTaylor(workspace.problem.objFun,guess.primal)
				# 	add_cuts!(guess,A=sparse(objApproximation.L'),loBs=[-Inf],upBs=[workspace.status.objUpB-objApproximation.c])
				# end
			end

			bestBound = solve!(guessList,workspace.nlpStepWS)
			if workspace.settings.verbose println("    Best Bound ",bestBound,".") end

			# update the list of feasible nodes found and update the objective upperbound
			feasibleIndices = findall(x->x.objUpB<workspace.settings.objectiveCutoff,guessList)
			if !isempty(feasibleIndices)
				workspace.status.objUpB = min(workspace.status.objUpB,bestBound)
				append!(workspace.feasibleNodes,guessList[feasibleIndices])
				collect_score_and_sort!(workspace.feasibleNodes,node->node.objUpB,algorithm=MergeSort)
			end
		end
		workspace.status.nlpTime += time() - timeCheckpoint


		################################ Termination 1 #########################
		# recompute optimality gaps
        (workspace.status.absoluteGap,workspace.status.relativeGap)=
			compute_gaps(workspace.status.objLoB,workspace.status.objUpB,workspace.settings.gapsComputationMode)

		# check optimality gaps satisfaction
		if workspace.status.absoluteGap < workspace.settings.absoluteGapTolerance ||
		   workspace.status.relativeGap < workspace.settings.relativeGapTolerance
		   if workspace.settings.verbose println(" HBB: Optimal Solution Found") end
			workspace.status.description = "optimalSolutionFound"
			break
		end

		# check if we have collected enough feasible points already
		if workspace.settings.numIncumbentsLimit > 0 && length(workspace.feasibleNodes) >= workspace.settings.numSolutionsLimit
			if workspace.settings.verbose println(" HBB: Interrupted") end
			workspace.status.description = "interrupted"
			break
		end



		######################### Linearizations Generation ########################
		timeCheckpoint = time()
		# generate new linear relaxations
		counter = get_numConstraints(workspace)
		newRelaxationsIDs = Int[]
		newRelaxationsSet = LinearConstraintSet(A=cnsJacobianType(zeros(0,numVars)),loBs=Float[],upBs=Float[])
		for (k,guess) in enumerate(guessList)
			if workspace.settings.limitedMemory
				cnssVals = evaluate(workspace.problem.cnsSet,guess.primal)
				activeCnssIndices = filter(i->max(cnssVals[i]-cnsUpBs[i],cnsLoBs[i]-cnssVals[i])>=-workspace.settings.primalTolerance,workspace.nonLinearCnssIndices)
				append!(newRelaxationsSet,linearRelaxation(workspace.problem.cnsSet,guess.primal)[activeCnssIndices])
				append!(newRelaxationsIDs,activeCnssIndices)
			else
				append!(newRelaxationsSet,linearRelaxation(workspace.problem.cnsSet,guess.primal)[workspace.nonLinearCnssIndices])
				append!(newRelaxationsIDs,workspace.nonLinearCnssIndices)
			end
		end
		workspace.status.rlxTime += time() - timeCheckpoint


		######################### MIP ########################
		timeCheckpoint = time()
		if workspace.settings.verbose
			println("- MILP step:")
		end

		# update the cutoff for the mip relaxation
		update_objectiveCutoff!(workspace.mipStepWS,get_suboptimality_threshold(workspace)-workspace.settings.primalTolerance)

		# insert new linearizations in the mip subproblem
		append_relaxations!(workspace.mipStepWS,newRelaxationsSet,newRelaxationsIDs)
		if workspace.settings.optimalControlInfo[1] > -1
			permutation = sortperm_oc(get_relaxationsSet(workspace.mipStepWS),workspace.settings.optimalControlInfo)
			permute_relaxations!(workspace.mipStepWS,permutation)
		end

		# solve the mip step and obtain new integer points and a new lower_bound
		mipStatus,MILPfeasibles,bestBound = solve!(workspace.mipStepWS)

		# create a list of guesses
		if !isempty(MILPfeasibles)

			# clean up the solutions
			for node in MILPfeasibles
				@. node.primal[dscIndices] = round(node.primal[dscIndices])
			end

			# remove the guesses that have been already considered
 			repeatedMask = repeatedDetector.(MILPfeasibles)

			# create the guess list by ignoring the repeated guesses
			guessList = MILPfeasibles[findall(iszero,repeatedMask)]

			# if the problem is mixed binary, insert the repeated guesses into the blacklist
			if is_mixedBinary(workspace.problem.varSet)
				repeatedInds = findall(!iszero,repeatedMask)
				for ind in repeatedInds
					addto_blackList!(workspace.mipStepWS,MILPfeasibles[ind].primal[dscIndices])
				end
			end
		else
			guessList = BBnode[]
		end


		if isempty(guessList)
			workspace.status.objLoB = workspace.status.objUpB
		else
			# update the lower bound
			workspace.status.objLoB = max(workspace.status.objLoB,bestBound)
			# consider only the first N(=num_processes) guesses (to avoid overloading the nlp step)
			deleteat!(guessList,workspace.settings.nlpProcesses+1:length(guessList))
		end

		if workspace.settings.verbose println("    # new Guesses found ",length(guessList),", Best Bound ",bestBound,".") end
		workspace.status.mipTime += time() - timeCheckpoint

		################################## Termination 2 ################################
		# recompute optimality gaps
		(workspace.status.absoluteGap,workspace.status.relativeGap)=
			compute_gaps(workspace.status.objLoB,workspace.status.objUpB,workspace.settings.gapsComputationMode)

		# check optimality gaps satisfaction
		if workspace.status.absoluteGap < workspace.settings.absoluteGapTolerance ||
		   workspace.status.relativeGap < workspace.settings.relativeGapTolerance
			if workspace.settings.verbose println(" HBB: Optimal Solution Found") end
			workspace.status.description = "optimalSolutionFound"
			break
		end

		# check for Infeasibilty
		if workspace.status.objLoB == Inf
			if workspace.settings.verbose println(" HBB: Infeasiblity Proven") end
			workspace.status.description = "infeasible"
			break
		end


		# check interruption rules
		if (workspace.settings.iterationsLimit>0 && numIterations == workspace.settings.iterationsLimit) ||
		   elapsedTime >= workspace.settings.timeLimit ||
		   workspace.settings.customInterruptionRule(workspace)
			if workspace.settings.verbose println(" HBB: Interrupted") end
			workspace.status.description = "interrupted"
			break
		end

		# if there are no new guesses and we haven't reached a terminal state: we are cycling
		# now, if the problem is mixed binary we can escape blacklisting the repeated assignments
		# in the mip step, otherwise, we need to throw an error
		if isempty(guessList) && !is_mixedBinary(workspace.problem.varSet)
			error("HBB: cycling")
		end

	end

	# store the time and the iterations required for the optimization
	workspace.status.numIterations += numIterations
	workspace.status.totalTime += time()-initialTime

	if workspace.settings.verbose
		println(" HBB termination:"*
	   		    " Obj.Bounds=("*string(workspace.status.objLoB)*","*string(workspace.status.objUpB)*")"*
	   			" Abs.Gap="*string(workspace.status.absoluteGap)*" Rel.Gap="*string(round(workspace.status.relativeGap,digits=2)))
    end
	return
end
