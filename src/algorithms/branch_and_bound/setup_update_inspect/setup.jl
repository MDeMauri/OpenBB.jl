# @Author: Massimo De Mauri <massimodemauri>
# @Date:   2019-02-06T16:22:41+01:00
# @Email:  massimo.demauri@gmail.com
# @Project: OpenBB
# @Filename: setup.jl
# @Last modified by:   massimo
# @Last modified time: 2021-02-09T22:21:58+01:00
# @License: LGPL-3.0
# @Copyright: {{copyright}}


function setup(problem::Problem,settings::BBsettings)::BBworkspace

	# load default settings
	if settings.subsolverSettings isa NullSubsolverSettings
		if problem.cnsSet isa LinearConstraintSet
			if problem.objFun isa LinearObjective
				settings.subsolverSettings = CLPsettings()
			elseif problem.objFun isa QuadraticObjective
				settings.subsolverSettings = CLPsettings()
			elseif problem.objFun isa ConvexObjective
				settings.subsolverSettings = IPOPTsettings()
			else
				error("Type of the problem not understood")
			end
		elseif problem.cnsSet isa ConvexConstraintSet
			if problem.objFun isa LinearObjective
				settings.subsolverSettings = IPOPTsettings()
			elseif problem.objFun isa QuadraticObjective
				settings.subsolverSettings = IPOPTsettings()
			elseif problem.objFun isa ConvexObjective
				settings.subsolverSettings = IPOPTsettings()
			else
				error("Type of the problem not understood")
			end
		else
			error("Type of the problem not understood")
		end
	end


	# check correctness of the setting and enforce their internal consistency
	enforce_consistency!(settings)

	if settings.numProcesses == 1 # single process_setup

		# construct the master BBworkspace
		problem_ = deepcopy(problem)
		settings_ = deepcopy(settings)
		workspace = setup_internal(problem_,settings_,NullSharedMemory(),drySetup=false)

	else # multi-process setup

		# add new processes if there aren't enough
		if nprocs() < settings.numProcesses
			addprocs(settings.numProcesses - nprocs())
		end

		# load OpenBB in the workers global scope
		@everywhere Main.eval(:(using OpenBB))
		workersList = workers()[1:settings.numProcesses-1]

		# synchronize the libraries with the remote workers if necessary
		libsManager.modified=false
		@sync for k in workersList
			@async remotecall_fetch(OpenBB.eval,k,:(remove_all_libs();load_libs($(list_libs()))))
		end

		# collect some data for the BBworkspace
		numVars = get_numVariables(problem)
		numCnss = get_numConstraints(problem)
		numDscVars = get_numDiscrete(problem)

		# construct the communication channels
		communicationChannels = Vector{BBnodeChannel}(undef,settings.numProcesses)
		for k in 1:settings.numProcesses
			communicationChannels[k] = BBnodeChannel(maxSerialSize_BBnode(numVars,numCnss,settings.maxNumberOfLocalCuts + is_mixedBinary(problem.varSet)))
            # obstruct the communication channels until the receiving worker is ready
			put!(communicationChannels[k],NullBBnode())
		end

		# construct shared Memory
		localObjLoBs = SharedVector{Float}(vcat([-Inf],repeat([Inf],settings.numProcesses-1)))
		globalObjUpB = SharedArray{Tuple{Float,Int8},1}([(Inf,0)])
		stats = SharedVector{Int}([0])
		arrestable = SharedVector{Bool}(repeat([false],settings.numProcesses))


		@sync if true # fake if for synchronization

			# build the remote workspaces
			for k in 1:length(workersList)
				@async remotecall_fetch(OpenBB.eval,workersList[k],
					:(workspace = OpenBB.setup_internal($problem,$settings,
												 	 OpenBB.BBsharedMemory($(communicationChannels[k]),
												 				  		   $(communicationChannels[k+1]),
																  		   $localObjLoBs,$globalObjUpB,
																  		   $stats,$arrestable),
												 	 drySetup=true)))
			end

			# construct the master BBworkspace
			problem_ = deepcopy(problem)
			settings_ = deepcopy(settings)
			workspace = setup_internal(problem_,settings_,
									BBsharedMemory(communicationChannels[end],
												   communicationChannels[1],
												   localObjLoBs,globalObjUpB,
												   stats,arrestable),
									drySetup=false)
		end

	end


    @sync return workspace

end


function setup_internal(problem::Problem,settings::BBsettings,sharedMemory::AbstractSharedMemory;drySetup::Bool=false)::BBworkspace

	# overwrite the subsolver settings depending on the branch and bound settings
	set_primalTolerance!(settings.subsolverSettings,min(settings.primalTolerance*1e-1,get_primalTolerance(settings.subsolverSettings)))
	set_dualTolerance!(settings.subsolverSettings,min(settings.dualTolerance,get_dualTolerance(settings.subsolverSettings)))
	 if settings.timeLimit < Inf
		 if get_timeLimit(settings.subsolverSettings) <= 0.0
			 set_timeLimit!(settings.subsolverSettings,settings.timeLimit)
		 else
			 set_timeLimit!(settings.subsolverSettings,min(settings.timeLimit,get_timeLimit(settings.subsolverSettings)))
		 end
	 end

	# actually create the workspace (with no root)
	subsolverWorkspace = setup(problem,settings.subsolverSettings,withPrecompilation=true,nodeType=BBnode)
	heuristicWorkspace = setup(problem,settings.heuristicsSettings,settings.subsolverSettings)
	workspace = BBworkspace(problem,subsolverWorkspace,sharedMemory,
							BBtree(),BBstatus(),BBupdatesRegister(),
							heuristicWorkspace,settings,
							BBblackList(get_numDiscrete(problem)),
							false,false,false)

	# initialize the pseudo costs
	initialize_pseudoCosts!(settings.pseudoCostsInitialization,workspace.problem.varSet.pseudoCosts)

	if drySetup
		# mark the workspace as empty
		workspace.status.description="empty"
		workspace.status.objLoB = Inf
	else
		# build the root node and insert it in the queue
		push_node!(workspace,:active,BBroot(workspace))
	end

	return workspace
end
