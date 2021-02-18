# @Author: Massimo De Mauri <massimo>
# @Date:   2019-05-23T11:28:34+02:00
# @Email:  massimo.demauri@gmail.com
# @Filename: BBjl
# @Last modified by:   massimo
# @Last modified time: 2021-01-25T17:53:19+01:00
# @License: LGPL-3.0
# @Copyright: {{copyright}}



mutable struct BBsettings <: AbstractSettings
    # settings for the subsolver
	subsolverSettings::SubsolverSettings					# those are passed directly to the subsolver
    # execution modifiers
    verbose::Bool                           				# print info during execution
    statusInfoPeriod::Float               					# frequency of status info print
    numProcesses::Int                       				# max number of processes to launch
    conservativismLevel::Int                  				# 0: normal BB, 1: suboptimal nodes are stored, 2: infeasible and suboptimal nodes are stored
	nodeSharingPeriod::Int									# send every n-th node to the neighbouring process
	# problem bounds
    primalTolerance::Float	                	    		# constraints violation tolerance
    dualTolerance::Float									# dual constraints violation tolerance
    objectiveCutoff::Float        	        				# look only for solutions that are better than the provided upper bound
	gapsComputationMode::Symbol     						# computes the optimality gaps from the objective bounds
    # priority rules
    expansionPriorityRule::Tuple{Symbol,Vararg}    			# ordering of the nodes in the activeQueue
	# optimality gaps computation
    branchingPriorityRule::Tuple{Symbol,Vararg}    			# ordering of the discrete variables for branching
    unreliablesPriority::Int      							# activeQueue insertion priority for unreliable nodes (-1->low, 0->normal, 1->high)
    # pseudo-costs
    pseudoCostsInitialization::Tuple{Symbol,Vararg}			# function returning the initialization of the pseudo-costs
    # stopping criteria
    customInterruptionRule::Function            			# user-defined stopping criterion
    timeLimit::Float                      					# max running time
    numSolutionsLimit::Int                  				# stop after the given number of integral solutions have been found
    absoluteGapTolerance::Float           					# stop if the absolute optimality gap is below the given level
    relativeGapTolerance::Float           					# stop if the relative optimality gap is below the given level
    # heuristic
	heuristicsSettings::BBheuristicsSettings				# defines which heuristics to use
    heuristicsTriggerCondition::Tuple{Symbol,Vararg}		# rule to trigger the heuristic
	# preprocessing
	withBoundsPropagation::Bool								# enables bounds propagation
	# algorithm modifiers
	optimalControlInfo::Tuple{Int,Int,Int}					# enables some optimal-control-specific rutines. the values are (numStates,numAlgebraicVars,numControls)
	acceptUnreliableSolutions::Bool							# consider solutions also in case of unreliability
	maxNumberOfLocalCuts::Int								# Each node carries around this number of cuts to improve the relaxations


	# constructor
	function BBsettings(;subsolverSettings::SubsolverSettings=NullSubsolverSettings(),
					     verbose::Bool=false,
	                     statusInfoPeriod::Float=1.0,
	                     numProcesses::Int=1,
	                     conservativismLevel::Int=0,
						 nodeSharingPeriod::Int=2,
	                     primalTolerance::Float=1e-4,
						 dualTolerance::Float=1e-6,
	                     objectiveCutoff::Float=Inf,
						 gapsComputationMode::Symbol=:onUpperBound,
	                     expansionPriorityRule::Tuple{Symbol,Vararg}=(:lower_pseudoObjective,),
	                     branchingPriorityRule::Tuple{Symbol,Vararg}=(:pseudoIncrements_geomean,2,1),
	                     unreliablesPriority::Int=0,
	                     pseudoCostsInitialization::Tuple{Symbol,Vararg}=(:initialize_to_constant!,1e-4),
	                     customInterruptionRule::Function=x->false,
	                     timeLimit::Float=Inf,
	                     numSolutionsLimit::Int=0,
	                     absoluteGapTolerance::Float=1e-3,
	                     relativeGapTolerance::Float=1e-3,
						 heuristicsSettings::BBheuristicsSettings=OpenBB.BBroundingHeuristicsSettings(),
						 heuristicsTriggerCondition::Tuple{Symbol,Vararg}=(:trigger_on_pseudoObjective,),
						 withBoundsPropagation::Bool=false,
						 optimalControlInfo::Tuple{Int,Int,Int}=(-1,-1,-1),
						 acceptUnreliableSolutions::Bool=false,
						 maxNumberOfLocalCuts::Int=5
	                     )::BBsettings


	    return new(	subsolverSettings,
				   	verbose,statusInfoPeriod,numProcesses,conservativismLevel,nodeSharingPeriod,
	            	primalTolerance,dualTolerance,objectiveCutoff,
				  	gapsComputationMode,expansionPriorityRule,
				  	branchingPriorityRule,unreliablesPriority,
                  	pseudoCostsInitialization,customInterruptionRule,
                  	timeLimit,numSolutionsLimit,absoluteGapTolerance,relativeGapTolerance,
				  	heuristicsSettings,heuristicsTriggerCondition,withBoundsPropagation,
				  	optimalControlInfo,acceptUnreliableSolutions,maxNumberOfLocalCuts)
	end
end


function enforce_consistency!(settings::BBsettings)::Nothing

	# check correctness of the inputs
	@assert settings.numProcesses >= 0
	@assert 0 <= settings.conservativismLevel <= 2
	@assert settings.primalTolerance >= 0.0 && settings.dualTolerance >= 0.0
	@assert settings.timeLimit > 0.0 && settings.numSolutionsLimit >= 0
	@assert settings.absoluteGapTolerance >= 0.0 && settings.relativeGapTolerance >= 0.0
	@assert settings.maxNumberOfLocalCuts > 0
	@assert settings.nodeSharingPeriod >= 2

	# auto select number of processes
	if settings.numProcesses == 0 settings.numProcesses = div(Sys.CPU_THREADS,2) end

	# check the correctness and viability of the optimal control info
	if any(settings.optimalControlInfo.>-1)
		@assert all(settings.optimalControlInfo.>=0) && any(settings.optimalControlInfo.>0)
		if !WITH_MPC_ADDON
			error("OpenBB: in order to handle optimal control info, the MPCaddon is necessary. MPCaddon not found.")
		end
	end

	# ensure consistency of tolerances
	set_primalTolerance!(settings.subsolverSettings,min(get_primalTolerance(settings.subsolverSettings),settings.primalTolerance))
	set_dualTolerance!(settings.subsolverSettings,min(get_dualTolerance(settings.subsolverSettings),settings.dualTolerance))

	return
end
