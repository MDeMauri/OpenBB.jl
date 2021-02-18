# @Author: Massimo De Mauri <massimodemauri>
# @Date:   2019-02-06T18:10:22+01:00
# @Email:  massimo.demauri@gmail.com
# @Project: OpenBB
# @Filename: solve!.jl
# @Last modified by:   massimo
# @Last modified time: 2021-01-28T16:41:33+01:00
# @License: LGPL-3.0
# @Copyright: {{copyright}}

# This is the main function called to solve a branch and bound problem
function solve!(workspace::BBworkspace)::Nothing

	if workspace.settings.verbose println("BB: Started...") end

	# collect info
	numProcesses = workspace.settings.numProcesses

	# update the workspace if necessary
	if workspace.outdated
		update!(workspace)
	end

	if numProcesses > 1

		# synchronize the possible external clibs
		if libsManager.modified
			libsManager.modified=false
			@sync for k in 2:numProcesses
				@async remotecall_fetch(OpenBB.eval,k,:(OpenBB.remove_all_libs();
													  OpenBB.load_libs($(list_libs()))))
			end
		end

		# synchronize the pseudo costs
		remotePseudoCosts = Vector{Tuple{Matrix{Float},Matrix{Int}}}(undef,numProcesses-1)
		@sync for (k,id) in enumerate(workers()[1:numProcesses-1])
			@async remotePseudoCosts[k] = remotecall_fetch(OpenBB.eval,id,:(workspace.problem.varSet.pseudoCosts))
		end
		for ps in remotePseudoCosts
			workspace.problem.varSet.pseudoCosts = (workspace.problem.varSet.pseudoCosts[1] + ps[1],
													workspace.problem.varSet.pseudoCosts[2] + ps[2])

		end
		workspace.problem.varSet.pseudoCosts = (workspace.problem.varSet.pseudoCosts[1]/numProcesses,
												workspace.problem.varSet.pseudoCosts[2])
		@sync for id in workers()[1:numProcesses-1]
			@async remotecall_fetch(OpenBB.eval,id,:(workspace.problem.varSet.pseudoCosts = $workspace.problem.varSet.pseudoCosts))
		end
	end

	@sync begin

		# start the remote branch and bound processes
		for k in workers()[1:numProcesses-1]
			@async remotecall_fetch(OpenBB.eval,k,:(try OpenBB.run!(workspace)
												    catch err
														  #collect error log
														  log = sprint(showerror,err,catch_backtrace())
														  push!(OpenBB.errorLogs,log)
														  # try to kill the other processes
														  OpenBB.put!(workspace.sharedMemory.outputChannel,KillerNode(0),mode=:forced)
														  error(log)
													end))
		end

		# start the local BB process
		try run!(workspace)
		catch err
			#collect error log
			log = sprint(showerror,err,catch_backtrace())
			push!(OpenBB.errorLogs,log)
			# try to kill the other processes
			if !(workspace.sharedMemory isa NullSharedMemory)
				put!(workspace.sharedMemory.outputChannel,KillerNode(0),mode=:forced)
			end
			# go on throwing the error
			error(log)
		end
	end


    ############################## termination ##############################
	# id of the current process
    processId = myid()
	# global status
	status = get_status(workspace)

	if status.somethingWrong
		status.description = "error"
		if workspace.settings.verbose && processId == 1
            print_status(workspace)
            println("BB: Forceful Stop Due to Error")
        end

	elseif status.absoluteGap < workspace.settings.absoluteGapTolerance ||
           status.relativeGap < workspace.settings.relativeGapTolerance

        status.description = "optimalSolutionFound"
        if workspace.settings.verbose && processId == 1
            print_status(workspace)
            println("BB: Optimal Solution Found")
        end

	elseif !all(@. getfield(getindex(workspace.tree.solutions,1),:reliable))

		status.description = "unreliableSolution"
		if workspace.settings.verbose && processId == 1
            print_status(workspace)
            println("BB: Unreliable Solutions Found")
        end

	elseif length(workspace.tree.active) == 0 && status.objUpB == Inf
		if status.cutoffActive
			status.description = "cutoffInfeasible"
			if workspace.settings.verbose && processId == 1
				print_status(workspace)
				println("BB: Infeasibilty Due to Cutoff")
			end
		elseif status.blackListActive
			status.description = "blackListInfeasible"
			if workspace.settings.verbose && processId == 1
				print_status(workspace)
				println("BB: Infeasibilty Due to Black-List")
			end
		else
			status.description = "infeasible"
			if workspace.settings.verbose && processId == 1
			  print_status(workspace)
			  println("BB: Infeasibilty Proven")
			end
		end

    else
        status.description = "interrupted"
        if workspace.settings.verbose && processId == 1
            print_status(workspace)
            println("BB: Interrupted")
        end
    end

	# propagate the reliability status and the status description,
	# finally obstruct the input channels to prevent further node sendings
	@sync if numProcesses > 1
		@async workspace.status.reliable = status.reliable
		@async workspace.status.description = status.description
		@async put!(workspace.sharedMemory.inputChannel,NullBBnode())
		for k in 2:numProcesses
			@async remotecall_fetch(OpenBB.eval,k,:(workspace.status.reliable = $status.reliable;
												  workspace.status.description = $status.description;
												  OpenBB.put!(workspace.sharedMemory.inputChannel,OpenBB.NullBBnode())))
		end
	else
		workspace.status.reliable = status.reliable
		workspace.status.description = status.description
	end
	return
end
