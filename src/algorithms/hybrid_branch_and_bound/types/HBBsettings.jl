# @Author: Massimo De Mauri <massimo>
# @Date:   2020-11-23T11:47:37+01:00
# @Email:  massimo.demauri@gmail.com
# @Filename: HBBjl
# @Last modified by:   massimo
# @Last modified time: 2021-02-22T15:17:39+01:00
# @License: LGPL-3.0
# @Copyright: {{copyright}}


mutable struct HBBsettings <: AbstractSettings
	# subsolvers settings
	nlpSettings::SubsolverSettings						# settings to be passed to the NLP subsolver
	mipSettings::BBsettings								# settings to be passed to the MIP subsolver
    # execution modifiers
    verbose::Bool                           			# print info during execution
    nlpProcesses::Int                       			# max number of processes for nlp
	mipProcesses::Int									# max number of processes for mip
	# definition of nlp and mip step
	nlpStepType::Tuple{Symbol,Vararg}					# defines the kind of nlp step to use
	mipStepType::Tuple{Symbol,Vararg}					# defines the kind of mip step to use
	# algorithm modifiers
	conservativismLevel::Int                  			# 0: normal BB, 1: suboptimal nodes are stored, 2: infeasible and suboptimal nodes are stored
    # problem bounds
	optimalControlInfo::Tuple{Int,Int,Int}				# enables some optimal-control-specific rutines. the values are (numStates,numAlgebraicVars,numControls)
	acceptUnreliableSolutions::Bool						# consider solutions also in case of unreliability
	limitedMemory::Bool 								# use all constraints for linearizations or only active ones
	# problem bounds
    primalTolerance::Float	                			# constraints violation tolerance
    dualTolerance::Float								# dual constraints violation tolerance
    objectiveCutoff::Float        	        			# look only for solutions that are better than the provided upper bound
	gapsComputationMode::Symbol							# computes the optimality gaps from the objective bounds
    # stopping criteria
    timeLimit::Float                      				# max running time
	iterationsLimit::Int								# max number of iterations
    numIncumbentsLimit::Int                  			# stop after the given number of integral solutions have been found
    absoluteGapTolerance::Float           				# stop if the absolute optimality gap is below the given level
    relativeGapTolerance::Float           				# stop if the relative optimality gap is below the given level
	customInterruptionRule::Function					# user-defined stopping criterion
	# stopping criterion for MIP step
	mipStepStoppingCriterion::Tuple{Symbol,Vararg}
	# heuristic
	heuristicsSettings::HBBheuristicsSettings 			# defines which heuristics to use
    heuristicsTriggerCondition::Tuple{Symbol,Vararg}    	# rule to trigger the heuristic


	# constructor
	function HBBsettings(;  # subsolvers settings
							nlpSettings::SubsolverSettings=IPOPTsettings(),
							mipSettings::BBsettings=BBsettings(),
							# execution modifiers
						    verbose::Bool=false,
						    nlpProcesses::Int=1,
							mipProcesses::Int=1,
							# definition of nlp and mip step
							nlpStepType::Tuple{Symbol,Vararg}=(:OAnlpStep,),
							mipStepType::Tuple{Symbol,Vararg}=(:OAmipStep,),
							# algorithm modifiers
							conservativismLevel::Int=0,
							optimalControlInfo::Tuple{Int,Int,Int}=(-1,-1,-1),
							acceptUnreliableSolutions::Bool=false,
							limitedMemory::Bool=false,
							# problem bounds
						    primalTolerance::Float=1e-4,
						    dualTolerance::Float=1e-6,
						    objectiveCutoff::Float=Inf,
							gapsComputationMode::Symbol=:onAverage,
						    # stopping criteria
						    timeLimit::Float=Inf,
							iterationsLimit::Int=0,
						    numIncumbentsLimit::Int=0,
						    absoluteGapTolerance::Float=1e-4,
						    relativeGapTolerance::Float=1e-6,
							customInterruptionRule::Function=workspace->false,
							# stopping criterion for MIP step
							mipStepStoppingCriterion::Tuple{Symbol,Vararg}=(:fullOptimality,),
						    # heuristic
							heuristicsSettings::HBBheuristicsSettings=NullHBBheuristicsSettings(),
						    heuristicsTriggerCondition::Tuple{Symbol,Vararg}=(:never_trigger,)
							)::HBBsettings


		return new(nlpSettings,mipSettings,verbose,nlpProcesses,mipProcesses,nlpStepType,mipStepType,
				   conservativismLevel,optimalControlInfo,acceptUnreliableSolutions,limitedMemory,
				   primalTolerance,dualTolerance,objectiveCutoff,gapsComputationMode,timeLimit,iterationsLimit,
				   numIncumbentsLimit,absoluteGapTolerance,relativeGapTolerance,customInterruptionRule,
				   mipStepStoppingCriterion,heuristicsSettings,heuristicsTriggerCondition)


	end
end

function enforce_consistency!(settings::HBBsettings)::Nothing
	# check correctness of the inputs
	@assert settings.nlpProcesses >= 0 && settings.mipProcesses >= 0
	@assert 0 <= settings.conservativismLevel <= 2
	@assert settings.primalTolerance >= 0.0 && settings.dualTolerance >= 0.0
	@assert settings.timeLimit > 0.0 && settings.numIncumbentsLimit >= 0
	@assert settings.absoluteGapTolerance >= 0.0 && settings.relativeGapTolerance >= 0.0


	# check the correctness and viability of the optimal control info
	if any(settings.optimalControlInfo.>-1)
		@assert all(settings.optimalControlInfo.>=0) && any(settings.optimalControlInfo.>0)
		if !WITH_MPC_ADDON
			error("OpenBB: in order to handle optimal control info, the MPCaddon is necessary. MPCaddon not found.")
		end
	end

	# load default settings for inner steps
	if settings.nlpSettings isa NullSettings
		settings.nlpSettings = IPOPTsettings()
	end

	# select number of processes if no indication is given by the user
	if settings.nlpProcesses == 0 settings.nlpProcesses = div(Sys.CPU_THREADS,2) end
	if settings.mipProcesses == 0 settings.mipProcesses = div(Sys.CPU_THREADS,2) end

	# ensure consistency in num. of processes
	settings.mipSettings.numProcesses = max(settings.mipSettings.numProcesses,settings.mipProcesses)

	# ensure consistency of tolerances
	set_primalTolerance!(settings.nlpSettings,min(get_primalTolerance(settings.nlpSettings),settings.primalTolerance))
	set_dualTolerance!(settings.nlpSettings,min(get_dualTolerance(settings.nlpSettings),settings.dualTolerance))
	settings.mipSettings.primalTolerance = min(settings.mipSettings.primalTolerance,settings.primalTolerance)
	settings.mipSettings.dualTolerance = min(settings.mipSettings.dualTolerance,settings.dualTolerance)
	settings.mipSettings.absoluteGapTolerance = settings.absoluteGapTolerance
	settings.mipSettings.relativeGapTolerance = settings.relativeGapTolerance
	settings.mipSettings.objectiveCutoff = settings.objectiveCutoff

	# enforce conservativismLevel >= 1 (HBB need) and ensure consistency in node handling
	settings.mipSettings.conservativismLevel = max(1,settings.conservativismLevel,settings.mipSettings.conservativismLevel)
	settings.mipSettings.optimalControlInfo = settings.optimalControlInfo
	settings.mipSettings.acceptUnreliableSolutions = settings.acceptUnreliableSolutions

	# decide when to stop the mip process
	if settings.mipStepStoppingCriterion[1] == :fullOptimality
		# do nothing
	elseif settings.mipStepStoppingCriterion[1] == :numOfGuesses
		settings.mipSettings.numSolutionsLimit = max(settings.mipSettings.numSolutionsLimit,settings.mipStepStoppingCriterion[2])
	elseif settings.mipStepStoppingCriterion[1] == :relativeGap
		settings.mipSettings.relativeGapTolerance = settings.mipStepStoppingCriterion[2]
	elseif settings.mipStepStoppingCriterion[1] == :absoluteGap
		settings.mipSettings.absoluteGapTolerance = settings.mipStepStoppingCriterion[2]
	else
		error("HBB: unknown MIP step stopping criterion")
	end

	return
end
