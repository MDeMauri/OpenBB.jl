# @Author: Massimo De Mauri <massimo>
# @Date:   2019-08-12T14:23:52+02:00
# @Email:  massimo.demauri@gmail.com
# @Filename: mpc_shift_utils.jl
# @Last modified by:   massimo
# @Last modified time: 2021-02-12T14:41:22+01:00
# @License: LGPL-3.0
# @Copyright: {{copyright}}

# performs a backward shift on the variables of a node
function shift_variables_backward!(node::BBnode,varShift::Int,bestSolution::Array{Float64,1},workspace::BBworkspace,
								   newUpdatedVars::BitArray{1},newUpdatedCnss::BitArray{1})::Bool

	# check if the best solution is feasible
	for k in 1:length(bestSolution)
		if bestSolution[k] < node.varLoBs[k] - workspace.settings.primalTolerance ||
		   bestSolution[k] > node.varUpBs[k] + workspace.settings.primalTolerance
			return false
		end
	end

	# collect info
	numStates = length(bestSolution) - varShift

    # shift backward the primal solution
	@. node.primal[1:numStates] = bestSolution[varShift+1:end]
    @. node.primal[numStates+1:end-varShift] = node.primal[varShift+numStates+1:end]

	# shift backward the dual solution
    @. node.bndDual[1:end-varShift] = node.bndDual[varShift+1:end]

	# impose the new initial state
	@. node.varLoBs[1:numStates] = bestSolution[varShift+1:end]
	@. node.varUpBs[1:numStates] = bestSolution[varShift+1:end]

    # shift backward the variable bounds
    @. node.varLoBs[numStates+1:end-varShift] = node.varLoBs[varShift+numStates+1:end]
    @. node.varUpBs[numStates+1:end-varShift] = node.varUpBs[varShift+numStates+1:end]

	# reset to default the bounds for the new variables
	variableBounds = get_bounds(workspace.problem.varSet)
	@. node.varLoBs[end-varShift+1:end] = variableBounds[1][end-varShift+1:end]
    @. node.varUpBs[end-varShift+1:end] = variableBounds[2][end-varShift+1:end]

	# shift the local cuts backward
	cutsOffset = node.cutSet.A[:,1:varShift]*bestSolution[1:varShift]
	node.cutSet.loBs -= cutsOffset
	node.cutSet.upBs -= cutsOffset
	node.cutSet.A = hcat(node.cutSet.A[:,varShift+1:end],spzeros(size(node.cutSet.A,1),varShift))

	# mark the initial state and the last variables for preprocessing
	newUpdatedVars[1:numStates] .= true
	newUpdatedVars[end-varShift+1:end] .= true

    return true
end

# performs a backward shift on the constraint dual of a node
function shift_constraints_backward!(node::BBnode,cnsShift::Int64,workspace::BBworkspace,
	newUpdatedVars::BitArray{1},newUpdatedCnss::BitArray{1})::Bool

    # shift backward the constraints bounds and Lagrangian multiplier
    @. node.cnsDual[1:end-cnsShift] = node.cnsDual[cnsShift+1:end]
    @. node.cnsLoBs[1:end-cnsShift] = node.cnsLoBs[cnsShift+1:end]
	@. node.cnsUpBs[1:end-cnsShift] = node.cnsUpBs[cnsShift+1:end]

	# reset to default the bounds of the new constraints
	cnsBounds = get_bounds(workspace.problem.cnsSet)
	@. node.cnsLoBs[end-cnsShift+1:end] = cnsBounds[1][end-cnsShift+1:end]
	@. node.cnsUpBs[end-cnsShift+1:end] = cnsBounds[2][end-cnsShift+1:end]

    return true
end
