# @Author: Massimo De Mauri <massimo>
# @Date:   2020-11-26T19:21:41+01:00
# @Email:  massimo.demauri@gmail.com
# @Filename: update_problem.jl
# @Last modified by:   massimo
# @Last modified time: 2021-02-16T13:41:42+01:00
# @License: LGPL-3.0
# @Copyright: {{copyright}}

# inserts constraints at the end of the constraint set
function append_constraints!(workspace::HBBworkspace,newCnsSet::ConstraintSet;
                             suppressErrors::Bool=false)::Nothing

    return insert_constraints!(workspace,newCnsSet,get_numConstraints(workspace)+1,suppressErrors=true)
end

# insert constraints at the specified point of the constraint set
function insert_constraints!(workspace::HBBworkspace,newCnsSet::ConstraintSet,insertionPoint::Int;
                             suppressErrors::Bool=false)::Nothing

	# check if it is possible to make changes
	if !suppressErrors && workspace.status.description != "new" && workspace.settings.conservativismLevel == 0
		error("HBB: In order to correctly perform a constraint insertion \"conservativismLevel\" should be set to at least 1")
	end

	# collect info
	numVars = get_size(workspace.problem.varSet)
	numCnss1 = get_size(workspace.problem.cnsSet)
	numCnss2 = get_size(newCnsSet)

    # insert constraints in the problem
    insert!(workspace.problem.cnsSet,newCnsSet,insertionPoint)

	# adapt constraints indices
	indicesMap = collect(1:numCnss1)
	for k in length(indicesMap):-1:1
		if indicesMap[k] >= insertionPoint
			indicesMap[k] += numCnss2
		else
			break
		end
	end

	# adapt the relaxations indices
	rename_relaxations!(workspace.mipStepWS,indicesMap)


	# insert the new indices
	newLinearIndices = findall(iszero,newCnsSet.hesSparsity)
	newNonLinearIndices = findall(!iszero,newCnsSet.hesSparsity)
	for k in length(workspace.linearCnssIndices):-1:1
		if workspace.linearCnssIndices[k] >= insertionPoint
			workspace.linearCnssIndices[k] += numCnss2
		else
			splice!(workspace.linearCnssIndices,k+1:0,newLinearIndices .+ (insertionPoint-1))
			break
		end
	end
	for k in length(workspace.nonLinearCnssIndices):-1:1
		if workspace.nonLinearCnssIndices[k] >= insertionPoint
			workspace.nonLinearCnssIndices[k] += numCnss2
		else
			splice!(workspace.nonLinearCnssIndices,k+1:0,newNonLinearIndices .+ (insertionPoint-1))
			break
		end
	end

	#TOREMOVE:
	@assert issorted(workspace.linearCnssIndices) && issorted(workspace.nonLinearCnssIndices)

    # create new linear relaxations
	if !isempty(newLinearIndices)
		newRelaxations = linearRelaxation(newCnsSet,zeros(numVars))[newLinearIndices]
		newRelaxationsIDs = copy(newLinearIndices.+(insertionPoint-1))
	else
		numVars = get_size(workspace.problem.varSet)
		if workspace.problem.cnsSet isa ConvexConstraintSet{<:AbstractMatrix,SpMatrix{Float}}
			newRelaxations = LinearConstraintSet(A=spzeros(0,numVars),loBs=Int[],upBs=[])
		else
			newRelaxations = LinearConstraintSet(A=zeros(0,numVars),loBs=Int[],upBs=[])
		end
		newRelaxationsIDs = Int[]
	end

	if !isempty(newNonLinearIndices) && !isempty(workspace.feasibleNodes)
		nlCnsSet = newCnsSet[newNonLinearIndices]
		for node in workspace.feasibleNodes
			cnsVal = evaluate(nlCnsSet,node.primal)
			active = findall(!iszero,cnsVal < nlCnsSet.loBs + workspace.settings.primalTolerance ||
					 				 cnsVal > nlCnsSet.loBs - workspace.settings.primalTolerance)
			if !isempty(active)
				append_constraints!(newRelaxations,linearRelaxation(nlCnsSet,node.primal)[active])
				append!(newRelaxationsIDs,newNonLinearIndices[active].+(insertionPoint-1))
			else
				#TODO: we might think about not removing the feasibleNodes that do not violate the new constraints
			end
		end
	end

	# insert the new linear relaxations in the problem
	append_relaxations!(workspace.mipStepWS,newRelaxations,newRelaxationsIDs)

    # declares the workspace and all subcomponents outdated
    make_outdated!(workspace)

	return
end


# removes the selected constraints
function remove_constraints!(workspace::HBBworkspace,indices::Union{Vector{Int},UnitRange{Int}};suppressErrors::Bool=false)::Nothing

	if !suppressErrors && workspace.status.description != "new"
		error("HBB: In orderd to correctly remove constraints, conservativism level should be set to 2. Please, re-solve the problem.")
	end

	# collect info
	IDs = collect(1:get_numConstraints(workspace))

	# ensure the sortedness of the indices
	indices = sort(indices)

	# change the problem definition
	remove_constraints!(workspace.problem.cnsSet,indices)

	# propagate the changes to the MIP stepd
	remove_constraints!(workspace.mipStepWS,indices,suppressErrors=suppressErrors)

	# update the relaxations map
	counter = 0
	start = 1
	for id in vcat(indices,[0])
		 for k in start:length(IDs)
			 if IDs[k] == id
				 IDs[k] = 0
				 start = k+1
				 counter += 1
				 break
			 else
				 IDs[k] -= counter
			 end
		 end
	end
	rename_relaxations!(workspace.mipStepWS,IDs)

	 # mark the workspace as outdated
	 make_outdated!(workspace)

	return
end

function permute_constraints!(workspace::HBBworkspace,permutation::Vector{Int};suppressErrors::Bool=false)::Nothing

	# apply changes to the problem
	permute!(workspace.problem.cnsSet,permutation)

	# apply changes to the nodes in the feasibleNodes
	for node in workspace.feasibleNodes
		# permute the constraints bounds and dual solution
		permute!(node.cnsDual,permutation)
		permute!(node.cnsLoBs,permutation)
		permute!(node.cnsUpBs,permutation)
	end

	# apply changes to the mip step
	permute_constraints!(workspace.mipStepWS,permutation)

	# propagate outdated state
	make_outdated!(workspace.nlpStepWS)
	make_outdated!(workspace.heuristicsWS)

    # yep no update for the MIP (nothing changes there)

	return
end

function update_bounds!(workspace::HBBworkspace;
						varLoBs::Vector{Float}=Vector{Float}(),
						varUpBs::Vector{Float}=Vector{Float}(),
                        cnsLoBs::Vector{Float}=Vector{Float}(),
                        cnsUpBs::Vector{Float}=Vector{Float}(),
                        suppressErrors::Bool=false)::Nothing

	# ensure the correctness of the input
	@assert length(varLoBs)==length(workspace.problem.varSet.loBs) || length(varLoBs)==0
	@assert length(varUpBs)==length(workspace.problem.varSet.upBs) || length(varUpBs)==0
	@assert length(cnsLoBs)==length(workspace.problem.cnsSet.loBs) || length(cnsLoBs)==0
	@assert length(cnsUpBs)==length(workspace.problem.cnsSet.upBs) || length(cnsUpBs)==0


	# check if it is possible to make changes
	if !suppressErrors && workspace.status.description != "new"
		# check if it is possible to perform a bounds modification
		if workspace.settings.conservativismLevel == 0
			error("HBB: In order to correctly perform a bound restriction \"conservativismLevel\" should be set to at least 1. Please, re-solve the problem")
		# check if a bounds relaxation was requested
		elseif (length(cnsLoBs) > 0 && any(@. cnsLoBs < workspace.problem.cnsSet.loBs)) ||
		   	   (length(cnsUpBs) > 0 && any(@. cnsUpBs > workspace.problem.cnsSet.upBs)) ||
		   	   (length(varLoBs) > 0 && any(@. varLoBs < workspace.problem.varSet.loBs)) ||
		   	   (length(varUpBs) > 0 && any(@. varUpBs > workspace.problem.varSet.upBs))
			error("HBB: in order to correctly relax the bounds of a problem, \"conservativismLevel\" should be set to 2. Please, re-solve the problem")
		end
	end



	# modify the problem definition
	update_bounds!(workspace.problem.varSet,loBs=varLoBs,upBs=varUpBs)
	update_bounds!(workspace.problem.cnsSet,loBs=cnsLoBs,upBs=cnsUpBs)

	# propagate the changes to the mip
	update_bounds_variables!(workspace.mipStepWS,loBs=varLoBs,upBs=varUpBs,suppressErrors=suppressErrors)
	if !(isempty(cnsLoBs)) && !(isempty(cnsUpBs))
		(oldCnsLoBs,oldCnsUpBs) = get_bounds(workspace.problem.cnsSet)
		update_bounds_constraints!(workspace.mipStepWS,oldLoBs=oldCnsLoBs,oldUpBs=oldCnsUpBs,newLoBs=cnsLoBs,newUpBs=cnsUpBs,suppressErrors=suppressErrors)
	elseif !(isempty(cnsLoBs))
		(oldCnsLoBs,~) = get_bounds(workspace.problem.cnsSet)
		update_bounds_constraints!(workspace.mipStepWS,oldLoBs=oldCnsLoBs,newLoBs=cnsLoBs,suppressErrors=suppressErrors)
	elseif !(isempty(cnsUpBs))
		(~,oldCnsUpBs) = get_bounds(workspace.problem.cnsSet)
		update_bounds_constraints!(workspace.mipStepWS,oldUpBs=oldCnsUpBs,newUpBs=cnsUpBs,suppressErrors=suppressErrors)
	end

	# make the workspace outdated
	make_outdated!(workspace)

	return
end


function integralize_variables!(workspace::HBBworkspace,indices::Union{Vector{Int},UnitRange{Int}};newSos1Groups::Vector{Int}=Int[],
                                suppressErrors::Bool=false,localOnly::Bool=false)::Nothing
	integralize!(workspace.problem.varSet,indices)
	integralize_variables!(workspace.mipStepWS,indices)
	make_outdated(workspace.heuristicsWS)
	return
end




#
function fix_variables!(workspace::HBBworkspace,indices::Vector{Int},values::Vector{Float};
                        removeFixedVariables::Bool=false,suppressErrors::Bool=false)::Nothing

	# check the correctness of the input
	@assert length(indices) == length(values)
	@assert all(@. workspace.problem.varSet.loBs[indices] <= values <= workspace.problem.varSet.upBs[indices])

	# fix variables in the problem
	fix_variables!(workspace.problem,indices,values,removeFixedVariables=removeFixedVariables)

    # propagate the changes to MIP step
	fix_variables!(workspace.mipStepWS,indices,removeFixedVariables=removeFixedVariables,suppressErrors=suppressErrors)

	# make workspace outdated
	make_outdated!(workspace)

    return
end





function append_problem!(workspace::HBBworkspace,newProblem::Problem;suppressErrors::Bool=false)::Nothing

	# check correctness of the input
	@assert newProblem.objFun isa QuadraticObjective || newProblem.objFun isa LinearObjective

	# make a copy of the problem to append
	newProblem_ = deepcopy(newProblem)

	# collect info
	numVars1 = get_numVariables(workspace.problem)
	numVars2 = get_numVariables(newProblem_)
	numCnss1 = get_numConstraints(workspace.problem)
	numCnss2 = get_numConstraints(newProblem_)

	# modify the variable set
	append!(workspace.problem.varSet,newProblem_.varSet)

	# modify the constraint set
	append_variables!(workspace.problem.cnsSet,numVars2)
	insert_variables!(newProblem_.cnsSet,numVars1,1)
	append!(workspace.problem.cnsSet,newProblem_.cnsSet)

	# modify the objective function
	append_variables!(workspace.problem.objFun,numVars2)
	insert_variables!(newProblem_.objFun,numVars1,1)
	add!(workspace.problem.objFun,newProblem_.objFun)

    # insert the new linear constraints in the constraint set of the MIP step
	newLinearIndices = findall(iszero,newProblem.cnsSet.hesSparsity)
	newRelaxations = linearRelaxation(problem.cnsSet,zeros(numVars2))[newLinearIndices]
	newRelaxedProblem = Problem(varSet=deepcopy(newProblem.varSet),
								cnsSet=newRelaxations,
								objFun=deepcopy(newProblem.objFun))
	append_problem!(workspace.mipStepWS,newRelaxedProblem,newLinearIndices.+numCnss1)

	# mark the workspace as outdatedF
	make_outdated!(workspace)

	return
end



function addto_blackList!(workspace::HBBworkspace,assignment::Vector{Float})::Nothing

	# check correctness of the input
	@assert length(assignment) == get_numDiscrete(workspace.problem)

	# clean up the assignment
	toInsert = round.(assignment)

	# make sure that the assignment is discrete enough
	for i in 1:length(assignment)
		if abs(assignment[i]-toInsert[i])>workspace.settings.primalTolerance
			error("HBB: Attempted the Insertion of a Non-Discrete Assignment into the Black-List")
		end
	end

	# insert in black-list
	insert!(workspace.blackList,toInsert)

	# possibly remove the assignment from the feasible nodes
	toRemove = Int[]
	for node in workspace.feasibleNodes
		if lookup(x,workspace.blackList)
			push!(toRemove,k)
		end
	end
	deleteat!(workspace.feasibleNodes,toRemove)

	# update the blacklist of the MIP step
	addto_blackList!(workspace.mipStepWS,assignment)

	return
end
