# @Author: Massimo De Mauri <massimo>
# @Date:   2019-02-19T12:14:43+01:00
# @Email:  massimo.demauri@gmail.com
# @Filename: update_problem.jl
# @Last modified by:   massimo
# @Last modified time: 2019-11-22T11:23:31+01:00
# @License: LGPL-3.0
# @Copyright: {{copyright}}

#
function append_constraints!(workspace::BBworkspace{T1,T2,T3},constraintSet::T4;
                             suppressWarnings::Bool=false,localOnly::Bool=false)::Nothing where T1<:Problem where T2<:AbstractWorkspace where T3<:AbstractSharedMemory where T4 <: AbstractConstraintSet

    return insert_constraints!(workspace,constraintSet,get_numConstraints(workspace)+1,
                               suppressWarnings=suppressWarnings,localOnly=localOnly)
end

#
function insert_constraints!(workspace::BBworkspace{T1,T2,T3},constraintSet::T4,index::Int;
                             suppressWarnings::Bool=false,localOnly::Bool=false)::Nothing where T1<:Problem where T2<:AbstractWorkspace where T3<:AbstractSharedMemory where T4 <: AbstractConstraintSet

    # call the same function on the other workers
    @sync if !localOnly && !(workspace.sharedMemory isa NullSharedMemory)
        # call the local version of the function on the remote workers
        for p in 2:workspace.settings.numProcesses
            @async remotecall_fetch(Main.eval,p,:(
				OpenBB.insert_constraints!(workspace,$constraintSet,$index,suppressWarnings=true,localOnly=true)
			))
        end

        # call the local version of the function on the current process
        insert_constraints!(workspace,constraintSet,index,suppressWarnings=suppressWarnings,localOnly=true)

    else

		# check if it is possible to make changes
		if !suppressWarnings && workspace.status.description != "new" && !workspace.settings.interactiveMode && myid() == 1
			@warn "In order to correctly manipulate the problem formulation, OpenBB must be run in interactive mode"
		end

		# change the problem definition
		insert!(workspace.problem.cnsSet,constraintSet,index)

		# mark the workspace as outdated
		make_outdated!(workspace)

		# propagate the changes to the nodes
		push!(workspace.updatesRegister,insert_constraints!,(index,deepcopy(get_bounds(constraintSet))))
    end

    return
end


#
function remove_constraints!(workspace::BBworkspace{T1,T2,T3},indices::Array{Int,1};
                             suppressWarnings::Bool=false,
                             localOnly::Bool=true)::Nothing where T1<:Problem where T2<:AbstractWorkspace where T3<:AbstractSharedMemory


    @sync if !localOnly && !(workspace.sharedMemory isa NullSharedMemory)
        # call the local version of the function on the remote workers
        for p in 2:workspace.settings.numProcesses
            @async remotecall_fetch(Main.eval,p,:(
				OpenBB.remove_constraints!(workspace,$indices,suppressWarnings=true,localOnly=true)
			))
        end

        # call the local version of the function on the main process
        remove_constraints!(workspace,indices,suppressWarnings=suppressWarnings,localOnly=true)

    else

		## check if it is possible to make changes
		if !suppressWarnings && workspace.status.description != "new"
			@warn "Removing constraints after some iterations is potentially destructive, We hope you know what you are doing."
		end

		# change the problem definition
		remove_constraints!(workspace.problem.cnsSet,indices)

		# mark the workspace as outdated
		make_outdated!(workspace)

        # propagate the changes to the nodes
		push!(workspace.updatesRegister,remove_constraints!,(copy(indices),))
    end

    return
end


#
function permute_constraints!(workspace::BBworkspace{T1,T2,T3},permutation::Array{Int,1};
                              suppressWarnings::Bool=false,localOnly::Bool=false)::Nothing where T1<:Problem where T2<:AbstractWorkspace where T3<:AbstractSharedMemory

    @sync if !localOnly && !(workspace.sharedMemory isa NullSharedMemory)
        # call the local version of the function on the remote workers
        for p in 2:workspace.settings.numProcesses
            @async remotecall_fetch(Main.eval,p,:(
				OpenBB.permute_constraints!(workspace,$permutation,suppressWarnings=true,localOnly=true)
			))
        end

        # call the local version of the function on the main process
        permute_constraints!(workspace,permutation,suppressWarnings=suppressWarnings,localOnly=true)

    else

		# modify the problem definition
		permute_constraints!(workspace.problem.cnsSet,permutation)

		# mark the subsolverWS as updated
		make_outdated!(workspace.subsolverWS)

        # propagate the changes to the nodes
		push!(workspace.updatesRegister,permute_constraints!,(copy(permutation),))

    end

    # mark the subsolver workspace as outdated
	make_outdated!(workspace.subsolverWS)

    return
end


#
function update_bounds!(workspace::BBworkspace{T1,T2,T3};
						varLoBs::Array{Float64,1}=Array{Float64,1}(),
						varUpBs::Array{Float64,1}=Array{Float64,1}(),
                        cnsLoBs::Array{Float64,1}=Array{Float64,1}(),
                        cnsUpBs::Array{Float64,1}=Array{Float64,1}(),
                        suppressWarnings::Bool=false,localOnly::Bool=false)::Nothing where T1<:Problem where T2<:AbstractWorkspace where T3<:AbstractSharedMemory


    # ensure the correctness of the input
	@assert length(varLoBs)==length(workspace.problem.varSet.loBs) || length(varLoBs)==0
	@assert length(varUpBs)==length(workspace.problem.varSet.upBs) || length(varUpBs)==0
	@assert length(cnsLoBs)==length(workspace.problem.cnsSet.loBs) || length(cnsLoBs)==0
	@assert length(cnsUpBs)==length(workspace.problem.cnsSet.upBs) || length(cnsUpBs)==0


    @sync if !localOnly && !(workspace.sharedMemory isa NullSharedMemory)
        # call the local version of the function on the remote workers
        for p in 2:workspace.settings.numProcesses
            @async remotecall_fetch(Main.eval,p,:(
				OpenBB.update_bounds!(workspace,cnsLoBs=$cnsLoBs,cnsUpBs=$cnsUpBs,varLoBs=$varLoBs,varUpBs=$varUpBs,
                                      suppressWarnings=true,localOnly=true)
			))
        end

        # call the local version of the function on the main process
        update_bounds!(workspace,cnsLoBs=cnsLoBs,cnsUpBs=cnsUpBs,varLoBs=varLoBs,varUpBs=varUpBs,
                       suppressWarnings=suppressWarnings,localOnly=true)

    else

		# check if it is possible to make changes
		if !suppressWarnings && workspace.status.description != "new"
			if !workspace.settings.interactiveMode
				@warn "In order to correctly manipulate the problem formulation, OpenBB must be run in interactive mode"
			# check if a bounds relaxation was requested
		elseif (length(cnsLoBs) > 0 && any(@. cnsLoBs < workspace.problem.cnsSet.loBs)) ||
			   (length(cnsUpBs) > 0 && any(@. cnsUpBs > workspace.problem.cnsSet.upBs)) ||
			   (length(varLoBs) > 0 && any(@. varLoBs < workspace.problem.varSet.loBs)) ||
			   (length(varUpBs) > 0 && any(@. varUpBs > workspace.problem.varSet.upBs))
				@warn "Relaxing the bounds after some iterations is potentially destructive, please make sure that the new bound set is more restrictive than the old one."
			end
		end

		# modify the problem definition
		update_bounds!(workspace.problem.varSet,loBs=varLoBs,upBs=varUpBs)
		update_bounds!(workspace.problem.cnsSet,loBs=cnsLoBs,upBs=cnsUpBs)

        # propagate the changes to the nodes
		push!(workspace.updatesRegister,update_bounds!,(copy(varLoBs),copy(varUpBs),copy(cnsLoBs),copy(cnsUpBs)))

		# mark the workspace as outdated
		make_outdated!(workspace)
    end

    return
end

#
function append_problem!(workspace::BBworkspace{T1,T2,T3},problem::Problem;suppressWarnings::Bool=false,localOnly::Bool=false)::Nothing where T1<:Problem where T2<:AbstractWorkspace where T3<:AbstractSharedMemory


    @sync if !localOnly && !(workspace.sharedMemory isa NullSharedMemory)
        # call the local version of the function on the remote workers
        for p in 2:workspace.settings.numProcesses
            @async remotecall_fetch(Main.eval,p,:(
				OpenBB.append_problem!(workspace,$problem,suppressWarnings=true,localOnly=true)
			))
        end
        # call the local version of the function on the main process
        append_problem!(workspace,problem,suppressWarnings=suppressWarnings,localOnly=true)

    else

		# check if it is possible to make changes
		if !suppressWarnings && workspace.status.description != "new" && !workspace.settings.interactiveMode && myid() == 1
			@warn "In order to correctly manipulate the problem formulation, OpenBB must be run in interactive mode"
		end

		# make a copy of the problem to append
		localProblem = deepcopy(problem)

		# collect info
		numVars1 = get_numVariables(workspace.problem)
		numVars2 = get_numVariables(localProblem)
		numCnss1 = get_numConstraints(workspace.problem)
		numCnss2 = get_numConstraints(localProblem)

		# modify the variable set
		append_variables!(workspace.problem.varSet,localProblem.varSet)

		# modify the constraint set
		append_variables!(workspace.problem.cnsSet,numVars2)
		insert_variables!(localProblem.cnsSet,numVars1,1)
		append!(workspace.problem.cnsSet,localProblem.cnsSet)

		# modify the objective function
		append_variables!(workspace.problem.objFun,numVars2)
		insert_variables!(localProblem.objFun,numVars1,1)
		add!(workspace.problem.objFun,localProblem.objFun)

		# mark the workspace as outdated
		make_outdated!(workspace)

        # propagate the changes to the nodes
		push!(workspace.updatesRegister,insert_variables!,(numVars1+1,deepcopy(problem.varSet.vals),deepcopy(get_bounds(problem.varSet))))
		push!(workspace.updatesRegister,insert_constraints!,(numCnss1+1,deepcopy(get_bounds(problem.cnsSet))))

    end

    return
end


#
function integralize_variables!(workspace::BBworkspace{T1,T2,T3},newDscIndices::Array{Int,1};newSos1Groups::Array{Int,1}=Int[],
                                suppressWarnings::Bool=false,localOnly::Bool=false)::Nothing where T1<:Problem where T2<:AbstractWorkspace where T3<:AbstractSharedMemory

	# check the correctness of the input
	@assert length(newSos1Groups) == 0 || length(newSos1Groups) == length(newDscIndices)

    @sync if !localOnly && !(workspace.sharedMemory isa NullSharedMemory)
        # call the local version of the function on the remote workers
        for p in 2:workspace.settings.numProcesses
            @async remotecall_fetch(Main.eval,p,:(
				OpenBB.integralize_variables!(workspace,$newDscIndices,newSos1Groups=$newSos1Groups,suppressWarnings=true,localOnly=true)
			))
        end

        # call the local version of the function on the main process
        integralize_variables!(workspace,newDscIndices,newSos1Groups=newSos1Groups,suppressWarnings=suppressWarnings,localOnly=true)

    else

        if !suppressWarnings && workspace.status.description != "new" && !workspace.settings.interactiveMode && myid() == 1
			# check if it is possible to make changes
			@warn "In order to correctly manipulate the problem formulation, OpenBB must be run in interactive mode"
	        # check correctness of the inputs
			intersection = intersect(newDscIndices,workspace.problem.varSet.dscIndices)
			if lenght(intersection) > 0
				@warn "Integralization of already integral variables requested."
			end
		end

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
		push!(workspace.updatesRegister,round_variable_bounds!,(copy(newDscIndices),))

    end

    return
end


#
function fix_variables!(workspace::BBworkspace{T1,T2,T3},indices::Array{Int,1};values::Array{Float64,1},
                                suppressWarnings::Bool=false,localOnly::Bool=false)::Nothing where T1<:Problem where T2<:AbstractWorkspace where T3<:AbstractSharedMemory

	# check the correctness of the input
	@assert length(newSos1Groups) == 0 || length(newSos1Groups) == length(newDscIndices)

    @sync if !localOnly && !(workspace.sharedMemory isa NullSharedMemory)
        # call the local version of the function on the remote workers
        for p in 2:workspace.settings.numProcesses
            @async remotecall_fetch(Main.eval,p,:(
				OpenBB.fix_variables!(workspace,$newDscIndices,newSos1Groups=$newSos1Groups,suppressWarnings=true,localOnly=true)
			))
        end

        # call the local version of the function on the main process
        fix_variables!(workspace,newDscIndices,newSos1Groups=newSos1Groups,suppressWarnings=suppressWarnings,localOnly=true)

    else

        if !suppressWarnings && workspace.status.description != "new" && !workspace.settings.interactiveMode && myid() == 1
			# check if it is possible to make changes
			@warn "In order to correctly manipulate the problem formulation, OpenBB must be run in interactive mode"
		end

        # propagate the changes to the nodes
		push!(workspace.updatesRegister,fix_variables!,(copy(indices),copy(values)))

    end

    return
end




# for debug... remove this stuff
function set_objective!(workspace::BBworkspace,newObjective::T;suppressWarnings::Bool=false,localOnly::Bool=false)::Nothing where T<:AbstractObjective
	workspace.problem.objFun.Q = newObjective.Q
	workspace.problem.objFun.L = newObjective.L
	make_outdated!(workspace)
	return
end
