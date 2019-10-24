# @Author: Massimo De Mauri <massimodemauri>
# @Date:   2019-02-06T18:10:22+01:00
# @Email:  massimo.demauri@gmail.com
# @Project: OpenBB
# @Filename: solve!.jl
# @Last modified by:   massimo
# @Last modified time: 2019-10-23T18:56:38+02:00
# @License: LGPL-3.0
# @Copyright: {{copyright}}


# include referred code
include("./branch_and_solve!.jl")
include("./insert_node.jl")
include("./run!.jl")


# This is the main function called to solve a branch and bound problem
function solve!(workspace::BBworkspace)::Nothing


	# update the workspace if necessary
	if workspace.outdated
		update!(workspace)
	end
	@sync if true
		# start the remote branch and bound processes
		if workspace.settings.numProcesses > 1
			for k in 2:workspace.settings.numProcesses
				@async remotecall_fetch(Main.eval,k,:(OpenBB.run!(workspace)))
			end
		end

		# start the local BB process
		run!(workspace)
 	end




    ############################## termination ##############################
	# id of the current process
    processId = myid()
	# global status
	status = get_status(workspace)

    if status.absoluteGap < workspace.settings.absoluteGapTolerance ||
       status.relativeGap < workspace.settings.relativeGapTolerance

        status.description = "optimalSolutionFound"
        if workspace.settings.verbose && processId == 1
            print_status(workspace)
            println(" Exit: Optimal Solution Found")
        end

	elseif status.cutoffActive

		status.description = "cutoffInfeasible"
        if workspace.settings.verbose && processId == 1
            print_status(workspace)
            println(" Exit: Infeasibilty Due to Cutoff")
        end


    elseif length(workspace.activeQueue) == 0 && status.objUpB == Inf

        status.description = "infeasible"
        if workspace.settings.verbose && processId == 1
            print_status(workspace)
            println(" Exit: Infeasibilty Detected")
        end

    else
        status.description = "interrupted"
        if workspace.settings.verbose && processId == 1
            print_status(workspace)
            println(" Exit: Interrupted")
        end
    end

	# propagate the reliability status
	workspace.status.reliable = status.reliable
	@sync if workspace.settings.numProcesses > 1
		for k in 2:workspace.settings.numProcesses
			@async remotecall_fetch(Main.eval,k,:(workspace.status.reliable = $status.reliable))
		end
	end

	# propagate the status description
	workspace.status.description = status.description
	@sync if workspace.settings.numProcesses > 1
		for k in 2:workspace.settings.numProcesses
			@async remotecall_fetch(Main.eval,k,:(workspace.status.description = $status.description))
		end
	end

	return
end
