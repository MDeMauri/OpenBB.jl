# @Author: Massimo De Mauri <massimo>
# @Date:   2020-04-07T11:57:16+02:00
# @Email:  massimo.demauri@gmail.com
# @Filename: mpc_problem_manipulation.jl
# @Last modified by:   massimo
# @Last modified time: 2021-02-12T14:45:23+01:00
# @License: LGPL-3.0
# @Copyright: {{copyright}}

function shift_pseudocosts_backward!(varSet::VariableSet,varShift::Int)::Nothing
	varSet.pseudoCosts[1][1:end-varShift,:] =  varSet.pseudoCosts[1][varShift+1:end,:]
	varSet.pseudoCosts[2][1:end-varShift,:] =  varSet.pseudoCosts[2][varShift+1:end,:]
	varSet.pseudoCosts[2][end-varShift+1:end,:] .= div.(varSet.pseudoCosts[2][end-varShift+1:end,:],3)
	return
end

function shift_backward!(cnsSet::LinearConstraintSet,varShift::Int,referenceSolution::Vector{Float})::Nothing

	# check input
	@assert length(referenceSolution) >= varShift

	# compute the effect of fixing the variables
	deltaBounds = view(cnsSet.A,:,1:varShift)*view(referenceSolution,1:varShift)
	cnsSet.loBs -= deltaBounds
	cnsSet.upBs -= deltaBounds

	# shift the constraint matrix
	if cnsSet.A isa SpMatrix
		cnsSet.A = hcat(cnsSet.A[:,varShift+1:end],spzeros(size(cnsSet.A,1),varShift))
	else
		cnsSet.A = hcat(cnsSet.A[:,varShift+1:end],zeros(size(cnsSet.A,1),varShift))
	end

	return
end


function shift_backward!(objFun::LinearObjective,varShift::Int,referenceSolution::Vector{Float})::Nothing

	# check input
	@assert length(referenceSolution) >= varShift

	# compute the effect of fixing the head variables
	objFun.c += view(objFun.L,1:varShift)'*view(referenceSolution,1:varShift)

	if objFun.L isa SpVector
		objFun.L = vcat(objFun.L[varShift+1:end],spzeros(varShift))
	else
		objFun.L = vcat(objFun.L[varShift+1:end],zeros(varShift))
	end

	return
end

function shift_backward!(objFun::QuadraticObjective,varShift::Int,referenceSolution::Vector{Float})::Nothing

	# collect info
	numVars = length(objFun.L)

	# compute the effect of fixing the head variables
	objFun.c += .5*referenceSolution[1:varShift]'*objFun.Q[1:varShift,1:varShift]*referenceSolution[1:varShift] + objFun.L[1:varShift]'*referenceSolution[1:varShift]

	# linear part
	if objFun.L isa SpVector
		objFun.L[1:end-varShift] = objFun.L[varShift+1:end] + view(objFun.Q,varShift+1:numVars,1:varShift)*sparse(view(referenceSolution,1:varShift))
		objFun.L[end-varShift+1:end] = spzeros(varShift)
	else
		objFun.L[1:end-varShift] = objFun.L[varShift+1:end] + view(objFun.Q,varShift+1:numVars,1:varShift)*view(referenceSolution,1:varShift)
		objFun.L[end-varShift+1:end] = zeros(varShift)
	end

	# quadratic part
	if objFun.Q isa SpMatrix
		objFun.Q[1:end-varShift,1:end-varShift] = objFun.Q[varShift+1:end,varShift+1:end]
		objFun.Q[:,end-varShift+1:end] = spzeros(numVars,varShift)
		objFun.Q[end-varShift+1:end,:] = spzeros(varShift,numVars)
	else
		objFun.Q[1:end-varShift,1:end-varShift] = objFun.Q[varShift+1:end,varShift+1:end]
		objFun.Q[:,end-varShift+1:end] = zeros(numVars,varShift)
		objFun.Q[end-varShift+1:end,:] = zeros(varShift,numVars)
	end

	return
end


function remove_terms!(objFun::LinearObjective,indices::Array{Int,1})::Nothing
	objFun.L[indices] .= 0.0
	return
end

function remove_terms!(objFun::QuadraticObjective,indices::Array{Int,1})::Nothing
	objFun.Q[indices,indices] .= 0.0
	objFun.L[indices] .= 0.0
	return
end
