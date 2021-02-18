# @Author: Massimo De Mauri <massimo>
# @Date:   2020-11-26T19:02:28+01:00
# @Email:  massimo.demauri@gmail.com
# @Filename: update_workspace.jl
# @Last modified by:   massimo
# @Last modified time: 2021-02-12T14:49:41+01:00
# @License: LGPL-3.0
# @Copyright: {{copyright}}



# it marks the workspace as outdated
function make_outdated!(workspace::HBBworkspace)::Nothing
	workspace.outdated = true
	make_outdated!(workspace.nlpStepWS)
	make_outdated!(workspace.mipStepWS)
	make_outdated!(workspace.heuristicsWS)
	return
end


# performs the required updates
function update!(workspace::HBBworkspace)::Nothing

	# adapt the status to the changes
	workspace.status.objUpB = Inf
	workspace.status.objLoB = -Inf
	workspace.status.absoluteGap = Inf
	workspace.status.relativeGap = Inf
	workspace.status.numIncumbents = 0
	workspace.status.reliable = true
	workspace.status.cutoffActive = false
	workspace.status.description = "interrupted"

	# empty the list of solutions
	empty!(workspace.feasibleNodes)

	# reset the black-list
	workspace.blackList = BBblackList(get_numDiscrete(workspace.problem))
	clear_blackList!(workspace.mipStepWS)

	# propagate the update
	if workspace.nlpStepWS.outdated update!(workspace.nlpStepWS) end
	if workspace.mipStepWS.outdated update!(workspace.mipStepWS) end
	if workspace.heuristicsWS.outdated update!(heuristicsWS) end

	return
end


# return the workspace to the initial state
function reset!(workspace::HBBworkspace)::Nothing

	# reset the state
	workspace.status = HBBstatus()

	# empty the list of solutions
	empty!(workspace.feasibleNodes)

	# generate a new nlp step workspace
	nlpStepWS = setup(workspace.problem,OpenBB.eval(Symbol(workspace.settings.nlpStepType[1],:Settings))(workspace.settings.nlpSettings,workspace.settings.nlpStepType[2:end]...),workspace.settings.nlpProcesses)

	# generate a new mip step workspace
	if !isempty(workspace.linearCnssIndices)
		relaxationsSet = linearRelaxation(workspace.problem.cnsSet,zeros(get_size(workspace.problem.varSet)))[workspace.linearCnssIndices]
		relaxationsIDs = workspace.linearCnssIndices
	else
		if problem.cnsSet isa ConvexConstraintSet{SpMatrix{Float}}
			relaxationsSet = LinearConstraintSet(A=SpMatrix{Float}(I,0,get_size(workspace.problem.varSet)),loBs=Float[],upBs=Float[])
		else
			relaxationsSet = LinearConstraintSet(A=Matrix{Float}(undef,0,get_size(workspace.problem.varSet)),loBs=Float[],upBs=Float[])
		end
		relaxationsIDs = Int[]
	end
	relaxation = Problem(varSet=deepcopy(workspace.problem.varSet),cnsSet=relaxationsSet,objFun=deepcopy(workspace.problem.objFun))
	mipStepWS = setup(relaxation,relaxationsIDs,OpenBB.eval(Symbol(workspace.settings.mipStepType[1],:Settings))(workspace.settings.mipSettings,workspace.settings.mipStepType[2:end]...))

	# generate a new workspace for the heuristics
	heuristicsWS = setup(workspace.problem,workspace.settings.heuristicsSettings,workspace.settings.nlpSettings)

	# reset the black-list
	workspace.blackList = BBblackList(get_numDiscrete(workspace))

	return
end



function Base.filter!(selector::Function,workspace::HBBworkspace,suppressErrors::Bool)::Nothing
	filter!(selector,workspace.mipStepWS,suppressErrors=suppressErrors)
	filter!(selector,workspace.feasibleNodes)
	workspace.status.numIncumbents = length(workspace.feasibleNodes)

	if workspace.status.numIncumbents
		workspace.status.objUpB = workspace.feasibleNodes[1]
	else
		workspace.status.objUpB = Inf
	end
	workspace.status.description = "interrupted"
	return
end



# this function inserts the given set of new constraints into the OpenBB workspace
# preserving the sparsity pattern typical of optimal control problems
function insert_constraints_oc!(workspace::HBBworkspace,constraintSet::ConstraintSet;suppressErrors::Bool=false)::Nothing

    # preliminary check
	@assert all(workspace.settings.optimalControlInfo.>-1)

	# collect info
	numVarsPerStep = sum(workspace.settings.optimalControlInfo)

    # append the constraints in the workspace
    append_constraints!(workspace,constraintSet)

    # sort the constraints in the HBB workspace and the relaxations
    sort_constraints_oc!(workspace)

    return
end



## Optimal Control Specific Functions

# this function sorts the constraints in the workspace in order to ensure the
# sparsity pattern typical of optimal control problems
function sort_constraints_oc!(workspace::HBBworkspace;suppressErrors::Bool=false)::Nothing

	# preliminary check
	@assert all(workspace.settings.optimalControlInfo.>-1)

	# collect info
	numVarsPerStep = sum(workspace.settings.optimalControlInfo)

	# sort the constraints according to the shooting structure
	if !issorted_oc(workspace.problem.cnsSet,workspace.settings.optimalControlInfo)
		permutation = sortperm_oc(workspace.problem.cnsSet,workspace.settings.optimalControlInfo)
		permute_constraints!(workspace,permutation)
	end

	# sort the linearizations in the mip step
	relaxationsSet = get_relaxationsSet(workspace.mipStepWS)
	if !issorted_oc(relaxationsSet,workspace.settings.optimalControlInfo)
		permutation = sortperm_oc(relaxationsSets,workspace.settings.optimalControlInfo)
		permute_relaxations!(workspace.mipStepWS,permutation)
	end

    return
end
