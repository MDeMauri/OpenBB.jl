# @Author: Massimo De Mauri <massimo>
# @Date:   2019-05-23T11:28:34+02:00
# @Email:  massimo.demauri@gmail.com
# @Filename: BBsettings.jl
# @Last modified by:   massimo
# @Last modified time: 2019-11-20T15:45:53+01:00
# @License: LGPL-3.0
# @Copyright: {{copyright}}


# define abstract and null types
abstract type AbstractSettings end; struct NullSettings  <: AbstractSettings  end


mutable struct BBsettings <: AbstractSettings
    # execution modifiers
    verbose::Bool                           # print info during execution
    statusInfoPeriod::Float64               # frequency of status info print
    numProcesses::Int                       # max number of processes to launch
    interactiveMode::Bool                   # store in memory suboptimal solutions and nodes to allow later updates
    # problem bounds
    primalTolerance::Float64                # constraint violation tolerance
    objectiveCutoff::Float64                # look only for solutions that are better than the provided upper bound
    # priority rules
    expansionPriorityRule::Tuple            # ordering of the nodes in the activeQueue
    branchingPriorityRule::Tuple            # ordering of the discrete variables for branching
    unreliablesPriority::Int      			# activeQueue insertion priority for unreliable nodes (-1->low, 0->normal, 1->high)
    # pseudo-costs
    pseudoCostsInitialization::Tuple        # function returning the initialization of the pseudo-costs
    # stopping criteria
    customStoppingRule::Function            # user-defined stopping criterion
    timeLimit::Float64                      # max running time
    numSolutionsLimit::Int                  # stop after the given number of integral solutions have been found
    absoluteGapTolerance::Float64           # stop if the absolute optimality gap is below the given level
    relativeGapTolerance::Float64           # stop if the relative optimality gap is below the given level
    # heuristic search
    roundingHeuristicsThreshold::Float64    # use the simple rounding heuristics whenever the average fractionality is under the threshold
	# preprocessing
	withBoundsPropagation::Bool
	# algorithm modifiers
	optimalControlMode::Bool				# enables some optimal-control-specific rutines
	acceptUnreliableSolutions::Bool			# consider solutions also in case of unreliability
	maxNumberOflocalCuts::Int				# Each node carries around this number of cuts to improve the relaxations
end



function BBsettings(;verbose::Bool=false,
                     statusInfoPeriod::Float64 = 1.,
                     numProcesses::Int=1,
                     interactiveMode::Bool=false,
                     primalTolerance::Float64=1e-4,
                     objectiveCutoff::Float64=Inf,
                     expansionPriorityRule::Tuple=(lower_pseudoObjective,),
                     branchingPriorityRule::Tuple=(pseudoIncrements_geomean,),
                     unreliablesPriority::Int=0,
                     pseudoCostsInitialization::Tuple=(initialize_to_constant!,1e20),
                     customStoppingRule::Function=x->false,
                     timeLimit::Float64=Inf,
                     numSolutionsLimit::Int=0,
                     absoluteGapTolerance::Float64=1e-4,
                     relativeGapTolerance::Float64=1e-6,
                     roundingHeuristicsThreshold::Float64 = -1.,
					 withBoundsPropagation::Bool=false,
					 optimalControlMode::Bool=false,
					 acceptUnreliableSolutions::Bool=false,
					 maxNumberOflocalCuts::Int=5
                     )::BBsettings


    # check correctness of the inputs
    @assert numProcesses>=0
	if numProcesses == 0
    	numProcesses = div(Sys.CPU_THREADS,2)
	end

    return BBsettings(verbose,statusInfoPeriod,numProcesses,interactiveMode,
                      primalTolerance,objectiveCutoff,
                      expansionPriorityRule,branchingPriorityRule,unreliablesPriority,
                      pseudoCostsInitialization,customStoppingRule,
                      timeLimit,numSolutionsLimit,absoluteGapTolerance,relativeGapTolerance,
					  roundingHeuristicsThreshold,withBoundsPropagation,
					  optimalControlMode,acceptUnreliableSolutions,maxNumberOflocalCuts)
end
