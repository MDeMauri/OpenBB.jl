# @Author: Massimo De Mauri <massimo>
# @Date:   2019-02-19T12:14:43+01:00
# @Email:  massimo.demauri@gmail.com
# @Filename: update_problem.jl
# @Last modified by:   massimo
# @Last modified time: 2021-02-16T13:41:24+01:00
# @License: LGPL-3.0
# @Copyright: {{copyright}}

# inserts constraints at the end of the constraint set
function append_constraints!(workspace::BBworkspace,constraintSet::T;
                             suppressErrors::Bool=false,localOnly::Bool=false)::Nothing where T <: ConstraintSet

    return insert_constraints!(workspace,constraintSet,get_numConstraints(workspace)+1,
                               suppressErrors=suppressErrors,localOnly=localOnly)
end

# insert constraints at the specified point of the constraint set
function insert_constraints!(workspace::BBworkspace,constraintSet::T,index::Int;
                             suppressErrors::Bool=false,localOnly::Bool=false)::Nothing where T <: ConstraintSet

    # call the same function on the other workers
    @sync if !localOnly && !(workspace.sharedMemory isa NullSharedMemory)
        # call the local version of the function on the remote workers
        for p in workers()[1:workspace.settings.numProcesses-1]
            @async remotecall_fetch(OpenBB.eval,p,:(
				OpenBB.insert_constraints!(workspace,$constraintSet,$index,suppressErrors=suppressErrors,localOnly=true)
			))
        end

        # call the local version of the function on the current process
        insert_constraints!(workspace,constraintSet,index,suppressErrors=suppressErrors,localOnly=true)

    else

		# check if it is possible to make changes
		if !suppressErrors && workspace.status.description != "new" && workspace.settings.conservativismLevel == 0 && myid() == 1
			error("BB: In order to correctly perform a constraint insertion \"conservativismLevel\" should be set to at least 1")
		end

		# change the problem definition
		insert!(workspace.problem.cnsSet,constraintSet,index)

		# mark the workspace as outdated
		make_outdated!(workspace,infeasiblesToRecover=false,invalidLowerBounds=false)

		# propagate the changes to the nodes
		push!(workspace.updatesRegister,insert_constraints!,(index,get_bounds(constraintSet)))
    end

    return
end


# removes the selected constraints
function remove_constraints!(workspace::BBworkspace,indices::Union{Vector{Int},UnitRange{Int}};
                             suppressErrors::Bool=false,
                             localOnly::Bool=false)::Nothing

    @sync if !localOnly && !(workspace.sharedMemory isa NullSharedMemory)
        # call the local version of the function on the remote workers
        for p in workers()[1:workspace.settings.numProcesses-1]
            @async remotecall_fetch(OpenBB.eval,p,:(
				OpenBB.remove_constraints!(workspace,$indices,suppressErrors=true,localOnly=true)
			))
        end

        # call the local version of the function on the main process
        remove_constraints!(workspace,indices,suppressErrors=suppressErrors,localOnly=true)

    else
		# limit case
		if isempty(indices) return end

		# check if it is possible to make changes
		if !suppressErrors && workspace.status.description != "new" && workspace.settings.conservativismLevel < 2
			error("BB: In order to correctly remove constraints after some iterations \"conservativismLevel\" should be set to 2")
		end

		# change the problem definition
		remove_constraints!(workspace.problem.cnsSet,indices)

		# mark the workspace as outdated
		make_outdated!(workspace,infeasiblesToRecover=true,invalidLowerBounds=true)

        # propagate the changes to the nodes
		push!(workspace.updatesRegister,remove_constraints!,(indices,))
    end

    return
end


#
function permute_constraints!(workspace::BBworkspace,permutation::Vector{Int};
                              suppressErrors::Bool=false,localOnly::Bool=false)::Nothing

    @sync if !localOnly && !(workspace.sharedMemory isa NullSharedMemory)
        # call the local version of the function on the remote workers
        for p in workers()[1:workspace.settings.numProcesses-1]
            @async remotecall_fetch(OpenBB.eval,p,:(
				OpenBB.permute_constraints!(workspace,$permutation,suppressErrors=true,localOnly=true)
			))
        end

        # call the local version of the function on the main process
        permute_constraints!(workspace,permutation,suppressErrors=suppressErrors,localOnly=true)

    else

		# modify the problem definition
		permute!(workspace.problem.cnsSet,permutation)

		# mark the subsolverWS and the heuristicsWS as updated
		make_outdated!(workspace.subsolverWS)
		make_outdated!(workspace.heuristicsWS)

        # propagate the changes to the nodes
		push!(workspace.updatesRegister,permute_constraints!,(permutation,))

    end
    return
end


#
function update_bounds!(workspace::BBworkspace;
						varLoBs::Vector{Float}=Vector{Float}(),
						varUpBs::Vector{Float}=Vector{Float}(),
                        cnsLoBs::Vector{Float}=Vector{Float}(),
                        cnsUpBs::Vector{Float}=Vector{Float}(),
                        suppressErrors::Bool=false,localOnly::Bool=false)::Nothing


    # ensure the correctness of the input
	@assert length(varLoBs)==length(workspace.problem.varSet.loBs) || length(varLoBs)==0
	@assert length(varUpBs)==length(workspace.problem.varSet.upBs) || length(varUpBs)==0
	@assert length(cnsLoBs)==length(workspace.problem.cnsSet.loBs) || length(cnsLoBs)==0
	@assert length(cnsUpBs)==length(workspace.problem.cnsSet.upBs) || length(cnsUpBs)==0


    @sync if !localOnly && !(workspace.sharedMemory isa NullSharedMemory)
        # call the local version of the function on the remote workers
        for p in workers()[1:workspace.settings.numProcesses-1]
            @async remotecall_fetch(OpenBB.eval,p,:(
				OpenBB.update_bounds!(workspace,cnsLoBs=$cnsLoBs,cnsUpBs=$cnsUpBs,varLoBs=$varLoBs,varUpBs=$varUpBs,
                                      suppressErrors=true,localOnly=true)
			))
        end

        # call the local version of the function on the main process
        update_bounds!(workspace,cnsLoBs=cnsLoBs,cnsUpBs=cnsUpBs,varLoBs=varLoBs,varUpBs=varUpBs,
                       suppressErrors=suppressErrors,localOnly=true)

    else

		# check if it is possible to make changes
		relaxedConstraints = false
		if !suppressErrors && workspace.status.description != "new"
			# check if a bounds relaxation was requested
			if (length(cnsLoBs) > 0 && any(@. cnsLoBs < workspace.problem.cnsSet.loBs)) ||
			   (length(cnsUpBs) > 0 && any(@. cnsUpBs > workspace.problem.cnsSet.upBs)) ||
			   (length(varLoBs) > 0 && any(@. varLoBs < workspace.problem.varSet.loBs)) ||
			   (length(varUpBs) > 0 && any(@. varUpBs > workspace.problem.varSet.upBs))
			    error("BB: In order to correctly relax the variable/constraint bounds after some iterations \"conservativismLevel\" should be set to 2")
				relaxedConstraints = true
			elseif workspace.settings.conservativismLevel == 0
				error("BB: In order to correctly restrict the variable/constraint bounds \"conservativismLevel\" should be set to at least 1")
			end
		end

		# modify the problem definition
		update_bounds!(workspace.problem.varSet,loBs=varLoBs,upBs=varUpBs)
		update_bounds!(workspace.problem.cnsSet,loBs=cnsLoBs,upBs=cnsUpBs)

        # propagate the changes to the nodes
		push!(workspace.updatesRegister,update_bounds!,(varLoBs,varUpBs,cnsLoBs,cnsUpBs))

		# mark the workspace as outdated
		make_outdated!(workspace,infeasiblesToRecover=relaxedConstraints,
								 invalidLowerBounds=relaxedConstraints)
    end

    return
end

#
function append_problem!(workspace::BBworkspace,problem::Problem;suppressErrors::Bool=false,localOnly::Bool=false)::Nothing


    @sync if !localOnly && !(workspace.sharedMemory isa NullSharedMemory)
        # call the local version of the function on the remote workers
        for p in workers()[1:workspace.settings.numProcesses-1]
            @async remotecall_fetch(OpenBB.eval,p,:(
				OpenBB.append_problem!(workspace,$problem,suppressErrors=true,localOnly=true)
			))
        end
        # call the local version of the function on the main process
        append_problem!(workspace,problem,suppressErrors=suppressErrors,localOnly=true)

    else

		# check if it is possible to make changes
		if !suppressErrors && workspace.status.description != "new" && workspace.settings.conservativismLevel == 0 && myid() == 1
			@warn "BB: In order to correctly perform the requested operation \"conservativismLevel\" should be set to at least 1"
		end

		# make a copy of the problem to append
		localProblem = deepcopy(problem)

		# collect info
		numVars1 = get_numVariables(workspace.problem)
		numVars2 = get_numVariables(localProblem)
		numCnss1 = get_numConstraints(workspace.problem)
		numCnss2 = get_numConstraints(localProblem)

		# modify the variable set
		append!(workspace.problem.varSet,localProblem.varSet)

		# modify the constraint set
		append_variables!(workspace.problem.cnsSet,numVars2)
		insert_variables!(localProblem.cnsSet,numVars1,1)
		append!(workspace.problem.cnsSet,localProblem.cnsSet)

		# modify the objective function
		append_variables!(workspace.problem.objFun,numVars2)
		insert_variables!(localProblem.objFun,numVars1,1)
		add!(workspace.problem.objFun,localProblem.objFun)

		# mark the workspace as outdated
		make_outdated!(workspace,infeasiblesToRecover=false,invalidLowerBounds=false)

        # propagate the changes to the nodes
		push!(workspace.updatesRegister,insert_variables!,(numVars1+1,problem.varSet.vals,get_bounds(problem.varSet)))
		push!(workspace.updatesRegister,insert_constraints!,(numCnss1+1,get_bounds(problem.cnsSet)))

    end

    return
end


#
function integralize_variables!(workspace::BBworkspace,newDscIndices::Vector{Int};newSos1Groups::Vector{Int}=Int[],
                                suppressErrors::Bool=false,localOnly::Bool=false)::Nothing

	# check the correctness of the input
	@assert length(newSos1Groups) == 0 || length(newSos1Groups) == length(newDscIndices)

    @sync if !localOnly && !(workspace.sharedMemory isa NullSharedMemory)
        # call the local version of the function on the remote workers
        for p in workers()[1:workspace.settings.numProcesses-1]
            @async remotecall_fetch(OpenBB.eval,p,:(
				OpenBB.integralize_variables!(workspace,$newDscIndices,newSos1Groups=$newSos1Groups,suppressErrors=true,localOnly=true)
			))
        end

        # call the local version of the function on the main process
        integralize_variables!(workspace,newDscIndices,newSos1Groups=newSos1Groups,suppressErrors=suppressErrors,localOnly=true)

    else
		# ensure consistency of the inputs
		if length(newSos1Groups) == 0
			newSos1Groups = repeat([0],length(newDscIndices))
		end

        # add the new discrete variables
        append!(workspace.problem.varSet.dscIndices,copy(newDscIndices))
        append!(workspace.problem.varSet.sos1Groups,copy(newSos1Groups))
        workspace.problem.varSet.pseudoCosts =
				(vcat(workspace.problem.varSet.pseudoCosts[1],repeat(sum(workspace.problem.varSet.pseudoCosts[1],dims=1)/size(workspace.problem.varSet.pseudoCosts[1],1),length(newDscIndices),1)),
                 vcat(workspace.problem.varSet.pseudoCosts[2],repeat([0 0],length(newDscIndices),1)))

        tmpPerm = sortperm(workspace.problem.varSet.dscIndices)
        permute!(workspace.problem.varSet.dscIndices,tmpPerm)
        permute!(workspace.problem.varSet.sos1Groups,tmpPerm)
        workspace.problem.varSet.pseudoCosts =
				(workspace.problem.varSet.pseudoCosts[1][tmpPerm,:],workspace.problem.varSet.pseudoCosts[2][tmpPerm,:])

        # propagate the changes to the nodes
		push!(workspace.updatesRegister,round_variable_bounds!,(newDscIndices,))

		# mark the workspace as outdated
		make_outdated!(workspace,infeasiblesToRecover=false,invalidLowerBounds=false)

    end

    return
end


#
function fix_variables!(workspace::BBworkspace,indices::Vector{Int},values::Vector{Float};
                        removeFixedVariables::Bool=false,suppressErrors::Bool=false,localOnly::Bool=false)::Nothing

	# check the correctness of the input
	@assert length(indices) == length(values)
	@assert all(@. workspace.problem.varSet.loBs[indices] <= values <= workspace.problem.varSet.upBs[indices])

    @sync if !localOnly && !(workspace.sharedMemory isa NullSharedMemory)
        # call the local version of the function on the remote workers
        for p in workers()[1:workspace.settings.numProcesses-1]
            @async remotecall_fetch(OpenBB.eval,p,:(
				OpenBB.fix_variables!(workspace,$newDscIndices,newSos1Groups=$newSos1Groups,
									  removeFixedVariables=$removeFixedVariables,suppressErrors=true,localOnly=true)
			))
        end

        # call the local version of the function on the main process
        fix_variables!(workspace,newDscIndices,newSos1Groups=newSos1Groups,
					  removeFixedVariables=removeFixedVariables,suppressErrors=suppressErrors,localOnly=true)

    else

        if !suppressErrors && workspace.status.description != "new" && workspace.settings.conservativismLevel == 0 && myid() == 1
			# check if it is possible to make changes
			error("BB: In order to correctly fix the variables value \"conservativismLevel\" should be set to at least 1")
		end

		if !all(-workspace.settings.primalTolerance + workspace.problem.varSet.loBs[indices] .<= values[indices] .<= workspace.problem.varSet.upBs[indices] + worspace.settings.primalTolerance)
			if !suppressErrors && workspace.status.description != "new" && workspace.settings.conservativismLevel < 2
			   	error("BB: In order to correctly fix the variables to a value outside of the previous bounds \"conservativismLevel\" should be set to 2")
			end
			boundsViolated = true
		else
			boundsViolated = false
		end


		# fix variables in the problem
		fix_variables!(workspace.problem,indices,values,removeFixedVariables=removeFixedVariables)

        # propagate the changes to the nodes
		push!(workspace.updatesRegister,fix_variables!,(indices,values),removeFixedVariables=removeFixedVariables)

		# mark the workspace as outdated
		make_outdated!(workspace,infeasiblesToRecover=boundsViolated,invalidLowerBounds=boundsViolated)

    end

    return
end




function update_objectiveCutoff!(workspace::BBworkspace,newCutoff::Float;
                                 suppressErrors::Bool=false,localOnly::Bool=false)::Nothing


	 # check the correctness of the input
     if !suppressErrors && workspace.status.cutoffActive && newCutoff > workspace.settings.objectiveCutoff
         @warn "BB: In order to correctly relax the cutoff after some iterations \"conservativismLevel\" should be set to at least 1"
     end

    @sync if !localOnly && !(workspace.sharedMemory isa NullSharedMemory)
        # call the local version of the function on the remote workers
        for p in workers()[1:workspace.settings.numProcesses-1]
            @async remotecall_fetch(OpenBB.eval,p,:(
                OpenBB.update_objectiveCutoff!(workspace,$newCutoff,suppressErrors=$suppressErrors,localOnly=true)
            ))
        end

        # call the function on the local worker
        update_objectiveCutoff!(workspace,newCutoff,suppressErrors=suppressErrors,localOnly=true)

    else

        # change the cutoff
        workspace.settings.objectiveCutoff = newCutoff

        if workspace.status.objLoB > newCutoff
        # declare the problem cutoff-infeasible and invalidate all solutions

			workspace.status.objUpB = workspace.status.absoluteGap = workspace.status.relativeGap = Inf
            workspace.status.objLoB = newCutoff
            workspace.status.cutoffActive = true
            workspace.status.description = "cutoffInfeasible"

            # move solutions into inactive
            if workspace.settings.conservativismLevel >= 1
                for (node,~) in workspace.tree.solutions
                    insert_node!(workspace,:suboptimals,node)
                end
            end
            empty!(workspace.tree.solutions)

            workspace.status.numSolutions = 0
            if !(workspace.sharedMemory isa NullSharedMemory)
                workspace.sharedMemory.stats[1] = 0
            end

        elseif workspace.status.objUpB > newCutoff
        # declare the problem interrupted and invalidate all solutions

            # update the status
            workspace.status.objUpB = workspace.status.absoluteGap = workspace.status.relativeGap = Inf
            workspace.status.cutoffActive = true
            workspace.status.description = "interrupted"

            # move solutions into inactive
            if workspace.settings.conservativismLevel >= 1
                for (node,~) in workspace.tree.solutions
                    insert_node!(workspace,:suboptimals,node)
                end
            end
            empty!(workspace.tree.solutions)

            workspace.status.numSolutions = 0
            if !(workspace.sharedMemory isa NullSharedMemory)
                workspace.sharedMemory.stats[1] = 0
            end
        else
        # invalidate the solutions that do not respect the new cutoff and keep the others
            solutionsToElim = Vector{Int}()
            for (k,(solution,~)) in enumerate(workspace.tree.solutions)
                if solution.objUpB > newCutoff - workspace.settings.primalTolerance
                    push!(solutionsToElim,k)
                    if workspace.settings.conservativismLevel >= 1
                        insert_node!(workspace,:suboptimals,solution)
                    end
                    workspace.status.numSolutions -= 1
                    if !(workspace.sharedMemory isa NullSharedMemory)
                        workspace.sharedMemory.stats[1] -= 1
                    end
                end
            end
            deleteat!(workspace.tree.solutions,solutionsToElim)
        end
    end

    return
end


function addto_blackList!(workspace::BBworkspace,assignment::Vector{Float};suppressErrors::Bool=false,localOnly::Bool=false)::Nothing

	# check correctness of the input
	@assert length(assignment) == get_numDiscrete(workspace.problem)
	if !is_mixedBinary(workspace.problem.varSet)
		error("BB: black-listing allowed only for mixed-binary problems")
	end

	@sync if !localOnly && !(workspace.sharedMemory isa NullSharedMemory)

		# call the local function in all the remote workers
		for p in workers()[1:workspace.settings.numProcesses-1]
            @async remotecall_fetch(OpenBB.eval,p,:(OpenBB.addto_blackList!(workspace,$assignment,suppressErrors=$suppressErrors,localOnly=true)))
        end

		# perform the operation in the current process
		addto_blackList!(workspace,assignment,suppressErrors=suppressErrors,localOnly=true)

	else
		# clean up the assignment
		toInsert = round.(assignment)

		# make sure that the assignment is discrete enough
		for i in 1:length(assignment)
			if abs(assignment[i]-toInsert[i])>workspace.settings.primalTolerance
				error("BB: Attempted the Insertion of a Non-Discrete Assignment into the Black-List")
			end
		end

		# insert in black-list
		insert!(workspace.blackList,toInsert)

		# remove the newly blacklisted solutions
		dscIndices = workspace.problem.varSet.dscIndices
		toRemove = Int[]
		newObjUpB = Inf
		for (k,(solution,~)) in enumerate(workspace.tree.solutions)
			if lookup(workspace.blackList,solution.primal[dscIndices])
				push!(toRemove,k)
			elseif newObjUpB > solution.objUpB
				newObjUpB = solution.objUpB
			end
		end
		append!(workspace.tree.blacklisted,workspace.tree.solutions[toRemove])
		deleteat!(workspace.tree.solutions,toRemove)
		workspace.status.objUpB = newObjUpB
	end

	return
end


function addto_blackList!(workspace::BBworkspace,assignments::Vector{Vector{Float}};suppressErrors::Bool=false,localOnly::Bool=false)::Nothing

	# check correctness of the input
	for ass in assignments
		@assert length(ass) == get_numDiscrete(workspace.problem)
	end
	if !is_mixedBinary(workspace.problem.varSet)
		error("BB: black-listing allowed only for mixed-binary problems")
	end

	@sync if !localOnly && !(workspace.sharedMemory isa NullSharedMemory)

		# call the local function in all the remote workers
		for p in workers()[1:workspace.settings.numProcesses-1]
            @async remotecall_fetch(OpenBB.eval,p,:(OpenBB.addto_blackList!(workspace,$assignments,suppressErrors=$suppressErrors,localOnly=true)))
        end

		# perform the operation in the current process
		addto_blackList!(workspace,assignments,suppressErrors=suppressErrors,localOnly=true)

	else
		toInsert = Vector{Vector{Float}}()
		for (k,a) in enumerate(assignments)
			# clean up assignment
			toInsert[k] = round.(a)
			# make sure that the assignment discrete enough
			for i in 1:length(a)
				if abs(a[i] - toInsert[k][i])>workspace.settings.primalTolerance
					error("BB: Attempted the Insertion of a Non-Discrete Assignment into the Black-List")
				end
			end
		end
		# insert in black-list
		insert!(workspace.blackList,toInsert)

		# remove the newly blacklisted solutions
		dscIndices = workspace.problem.varSet.dscIndices
		toRemove = Int[]
		newObjUpB = Inf
		for (k,(solution,~)) in enumerate(workspace.tree.solutions)
			if lookup(workspace.blackList,solution.primal[dscIndices])
				push!(toRemove,k)
			elseif newObjUpB > solution.objUpB
				newObjUpB = solution.objUpB
			end
		end
		append!(workspace.tree.blacklisted,workspace.tree.solutions[toRemove])
		deleteat!(workspace.tree.solutions,toRemove)
		workspace.status.objUpB = newObjUpB
	end

	return
end


function clear_blackList!(workspace::BBworkspace;suppressErrors::Bool=false,localOnly::Bool=false)::Nothing

	@sync if !localOnly && !(workspace.sharedMemory isa NullSharedMemory)

		# call the local function in all the remote workers
		for p in workers()[1:workspace.settings.numProcesses-1]
	        @async remotecall_fetch(OpenBB.eval,p,:(OpenBB.clear_blackList!(workspace,suppressErrors=$suppressErrors,localOnly=true)))
	    end

		# perform the operation in the current process
		clear_blackList!(workspace,suppressErrors=suppressErrors,localOnly=true)

	else

		# put all blacklisted nodes back in the tree
		for (node,~) in workspace.tree.blacklisted
			insert_node!(workspace,:active,node)
		end
		empty!(workspace.tree.blacklisted)

		# generate a new empty black-list
		workspace.blackList = BBblackList(get_numDiscrete(workspace.problem))

	end

	return
end


## Optimal Control Specific Functions
function sort_constraints_oc!(problem::Problem,optimalControlInfo::Tuple{Int,Int,Int})::Nothing
	permutation = sortperm_oc(problem.cnsSet,optimalControlInfo)
	if !issorted(permutation)
		permute!(problem.cnsSet,permutation)
	end
	return
end
