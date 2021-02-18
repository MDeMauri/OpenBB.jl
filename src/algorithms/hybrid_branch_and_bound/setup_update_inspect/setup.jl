# @Author: Massimo De Mauri <massimo>
# @Date:   2020-11-23T16:17:13+01:00
# @Email:  massimo.demauri@gmail.com
# @Filename: setup.jl
# @Last modified by:   massimo
# @Last modified time: 2021-01-28T16:33:38+01:00
# @License: LGPL-3.0
# @Copyright: {{copyright}}


function setup(problem::Problem,settings::HBBsettings)::HBBworkspace

	# check inputs
	@assert settings.nlpProcesses >= 0 && settings.mipProcesses >= 0
	if problem.cnsSet isa LinearConstraintSet ||
	   !(problem.objFun isa LinearObjective || problem.objFun isa QuadraticObjective)

		error("Hybrid Branch and Bound can be used only for problems having\n
			   a non-linear constraint set and a linear/quadratic objective\n
			   use non-linear Branch and Bound or reformulate the problem")
	end

	# check correctness of settings and enforce their internal consistency
	enforce_consistency!(settings)

	# prepare multi-process setup
	numProcesses = max(settings.nlpProcesses,settings.mipProcesses)
	if numProcesses > 1
		# add new processes if there aren't enough
		if nprocs() < numProcesses
			addprocs(numProcesses - nprocs())
		end

		# load OpenBB in the workers global scope
		@everywhere Main.eval(:(using OpenBB))

		# synchronize the libraries with the remote workers if necessary
		libsManager.modified=false
		@sync for k in workers()[1:numProcesses-1]
			@async remotecall_fetch(OpenBB.eval,k,:(OpenBB.remove_all_libs();OpenBB.load_libs($(list_libs()))))
		end
	end

	# store a local version of the input problem
	problem_ = Problem(varSet=deepcopy(problem.varSet),cnsSet=deepcopy(problem.cnsSet),objFun=deepcopy(problem.objFun))

	# generate the nlp step workspace
	nlpStepWS = setup(problem_,OpenBB.eval(Symbol(settings.nlpStepType[1],:Settings))(settings.nlpSettings,settings.nlpStepType[2:end]...),settings.nlpProcesses)

	# generate the mip step workspace
	linearCnssIndices = findall(iszero,problem.cnsSet.hesSparsity)
	nonLinearCnssIndices = collect(1:get_numConstraints(problem)); deleteat!(nonLinearCnssIndices,linearCnssIndices)
	if !isempty(linearCnssIndices)
		relaxationsSet = linearRelaxation(problem.cnsSet,zeros(get_numVariables(problem)))[linearCnssIndices]
		relaxationsIDs = copy(linearCnssIndices)
	else
		if problem.cnsSet isa ConvexConstraintSet{<:Any,SpMatrix{Float}}
			relaxationsSet = LinearConstraintSet(A=spzeros(0,get_size(problem.varSet)),loBs=Float[],upBs=Float[])
		else
			relaxationsSet = LinearConstraintSet(A=zeros(0,get_size(problem.varSet)),loBs=Float[],upBs=Float[])
		end
		relaxationsIDs = Int[]
	end
	relaxation = Problem(varSet=deepcopy(problem.varSet),cnsSet=relaxationsSet,objFun=deepcopy(problem.objFun))
	settings.mipSettings.conservativismLevel = max(settings.mipSettings.conservativismLevel,1)
	mipStepWS = setup(relaxation,relaxationsIDs,OpenBB.eval(Symbol(settings.mipStepType[1],:Settings))(settings.mipSettings,settings.mipStepType[2:end]...))

	# generate workspace for heuristics
	heuristicsWS = setup(problem_,settings.heuristicsSettings,settings.nlpSettings)

	return HBBworkspace(problem_,nlpStepWS,mipStepWS,BBnode[],linearCnssIndices,nonLinearCnssIndices,HBBstatus(),heuristicsWS,settings,
						BBblackList(get_numDiscrete(problem_)),false)
end
