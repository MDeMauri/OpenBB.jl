# @Author: Massimo De Mauri <massimo>
# @Date:   2020-10-22T13:33:52+02:00
# @Email:  massimo.demauri@gmail.com
# @Filename: IPOPT_interface.jl
# @Last modified by:   massimo
# @Last modified time: 2021-01-27T17:39:08+01:00
# @License: LGPL-3.0
# @Copyright: {{copyright}}

using Ipopt

## Settings ##########################################################
# structure used to hold the settings for IPOPT
mutable struct IPOPTsettings <: SubsolverSettings
	# generic options
	unlistedOptions::Dict{String,Union{String,Int,Float}}
	# convergence
	warm_start_init_point::String
	acceptable_compl_inf_tol::Float
	acceptable_constr_viol_tol::Float
	acceptable_dual_inf_tol::Float
	acceptable_iter::Int
	acceptable_obj_change_tol::Float
	acceptable_tol::Float
	compl_inf_tol::Float
	constr_viol_tol::Float
	diverging_iterates_tol::Float
	dual_inf_tol::Float
	tol::Float
	max_cpu_time::Float
	max_iter::Int
	mu_target::Float
	s_max::Float
	# linear solver
	linear_solver::String
	linear_system_scaling::String
	# NLP
	mu_oracle::String
	mu_strategy::String
	bound_relax_factor::Float
	check_derivatives_for_naninf::String
	dependency_detection_with_rhs::String
	dependency_detector::String
	fixed_variable_treatment::String
	hessian_constant::String
	honor_original_bounds::String
	jac_c_constant::String
	jac_d_constant::String
	kappa_d::Float
	nlp_lower_bound_inf::Float
	nlp_upper_bound_inf::Float
	num_linear_variables::Int
	# NLP scaling
	nlp_scaling_constr_target_gradient::Float
	nlp_scaling_max_gradient::Float
	nlp_scaling_method::String
	nlp_scaling_min_value::Float
	nlp_scaling_obj_target_gradient::Float
	obj_scaling_factor::Float
	# output
	print_frequency_iter::Int
	print_frequency_time::Float
	print_level::Int # the actual default is 5
	print_timing_statistics::String
	print_user_options::String
	# restoration phase
	bound_mult_reset_threshold::Float
	constr_mult_reset_threshold::Float
	evaluate_orig_obj_at_resto_trial::String
	expect_infeasible_problem::String
	expect_infeasible_problem_ctol::Float
	expect_infeasible_problem_ytol::Float
	max_resto_iter::Int
	max_soft_resto_iters::Int
	required_infeasibility_reduction::Float
	resto_failure_feasibility_threshold::Float
	resto_penalty_parameter::Float
	resto_proximity_weight::Float
	soft_resto_pderror_reduction_factor::Float
	start_with_resto::String
end

function IPOPTsettings(;
	# generic options
	unlistedOptions::Dict{String,Union{String,Int,Float}} = Dict{String,Union{String,Int,Float}}(),
	# convergence
	warm_start_init_point::String="yes",
	acceptable_compl_inf_tol::Float=0.01,
	acceptable_constr_viol_tol::Float=0.01,
	acceptable_dual_inf_tol::Float=1e10,
	acceptable_iter::Int=15,
	acceptable_obj_change_tol::Float=1e20,
	acceptable_tol::Float=1e-6,
	compl_inf_tol::Float=1e-4,
	constr_viol_tol::Float=1e-4,
	diverging_iterates_tol::Float=1e20,
	dual_inf_tol::Float=1e0,
	tol::Float=1e-8,
	max_cpu_time::Float=1e6,
	max_iter::Int=3000,
	mu_target::Float=0.0,
	s_max::Float=1e2,
	# linear solver
	linear_solver::String="ma27",
	linear_system_scaling::String="mc19",
	# NLP
	mu_oracle::String="probing", # actual default: "quality-function"
	mu_strategy::String="adaptive", # actual default: "monotone"
	bound_relax_factor::Float=0.0, # actual default: 1e-8
	check_derivatives_for_naninf::String="no",
	dependency_detection_with_rhs::String="no",
	dependency_detector::String="none",
	fixed_variable_treatment::String="make_parameter",
	hessian_constant::String="no",
	honor_original_bounds::String="yes",
	jac_c_constant::String="no",
	jac_d_constant::String="no",
	kappa_d::Float=1e-5,
	nlp_lower_bound_inf::Float=-1e19,
	nlp_upper_bound_inf::Float=1e19,
	num_linear_variables::Int=0,
	# NLP scaling
	nlp_scaling_constr_target_gradient::Float=0.0,
	nlp_scaling_max_gradient::Float=1e2,
	nlp_scaling_method::String="gradient-based",
	nlp_scaling_min_value::Float=1e-8,
	nlp_scaling_obj_target_gradient::Float=0.0,
	obj_scaling_factor::Float=1e0,
	# output
	print_frequency_iter::Int=1,
	print_frequency_time::Float=0.0,
	print_level::Int=0, # the actual default is 5
	print_timing_statistics::String="no",
	print_user_options::String="no",
	# restoration phase
	bound_mult_reset_threshold::Float=1e4,
	constr_mult_reset_threshold::Float=0.0,
	evaluate_orig_obj_at_resto_trial::String="yes",
	expect_infeasible_problem::String="no",
	expect_infeasible_problem_ctol::Float=1e-3,
	expect_infeasible_problem_ytol::Float=1e8,
	max_resto_iter::Int=3000000,
	max_soft_resto_iters::Int=10,
	required_infeasibility_reduction::Float=9e-1,
	resto_failure_feasibility_threshold::Float=0.0,
	resto_penalty_parameter::Float=1e3,
	resto_proximity_weight::Float=1e0,
	soft_resto_pderror_reduction_factor::Float=1.0-1e-4,
	start_with_resto::String="no"
	)

    return IPOPTsettings(unlistedOptions,warm_start_init_point,acceptable_compl_inf_tol,acceptable_constr_viol_tol,acceptable_dual_inf_tol,
						 acceptable_iter,acceptable_obj_change_tol,acceptable_tol,compl_inf_tol,constr_viol_tol,
						 diverging_iterates_tol,dual_inf_tol,tol,max_cpu_time,max_iter,mu_target,s_max,linear_solver,
						 linear_system_scaling,mu_oracle,mu_strategy,bound_relax_factor,check_derivatives_for_naninf,
						 dependency_detection_with_rhs,dependency_detector,fixed_variable_treatment,hessian_constant,
						 honor_original_bounds,jac_c_constant,jac_d_constant,kappa_d,nlp_lower_bound_inf,nlp_upper_bound_inf,
						 num_linear_variables,nlp_scaling_constr_target_gradient,nlp_scaling_max_gradient,nlp_scaling_method,
						 nlp_scaling_min_value,nlp_scaling_obj_target_gradient,obj_scaling_factor,print_frequency_iter,
						 print_frequency_time,print_level,print_timing_statistics,print_user_options,bound_mult_reset_threshold,
						 constr_mult_reset_threshold,evaluate_orig_obj_at_resto_trial,expect_infeasible_problem,
						 expect_infeasible_problem_ctol,expect_infeasible_problem_ytol,max_resto_iter,max_soft_resto_iters,
						 required_infeasibility_reduction,resto_failure_feasibility_threshold,resto_penalty_parameter,
						 resto_proximity_weight,soft_resto_pderror_reduction_factor,start_with_resto)
end


## Utility Functions ##########################################################

# wrapper to get verbosity of each solver in the same way
function get_verbosity(settings::IPOPTsettings)::Bool
	return settings.print_level > 0
end

# wrapper to get primal tolerance of each solver in the same way
function get_primalTolerance(settings::IPOPTsettings)::Float
	return settings.constr_viol_tol
end

# wrapper to get dual tolerance of each solver in the same way
function get_dualTolerance(settings::IPOPTsettings)::Float
	return settings.dual_inf_tol
end

# wrapper to get iterations limit of each solver in the same way
function get_iterationsLimit(settings::IPOPTsettings)::Int
	return settings.max_iter
end

# wrapper to get time limit of each solver in the same way
function get_timeLimit(settings::IPOPTsettings)::Float
	return settings.max_cpu_time
end

# wrapper to set verbosity of each solver in a consistent way
function set_verbosity!(settings::IPOPTsettings,verbose::Bool)::Nothing
	if verbose
		settings.print_level = 3
	end
	return
end

# wrappers to set primal tolerance of each solve in a consistent way
function set_primalTolerance!(settings::IPOPTsettings,tolerance::Float)::Nothing
	settings.constr_viol_tol = tolerance
	return
end

# wrappers to set dual tolerance of each solve in a consistent way
function set_dualTolerance!(settings::IPOPTsettings,tolerance::Float)::Nothing
	settings.dual_inf_tol = tolerance
	return
end

# wrappers to set the iteration limit of each solver in a consistent way
function set_iterationsLimit!(settings::IPOPTsettings,limit::Int)::Nothing
	settings.max_iter=limit
	return
end

# wrappers to set the time limit of each solver in a consistent way
function set_timeLimit!(settings::IPOPTsettings,limit::Float)::Nothing
	settings.max_cpu_time=limit
	return
end

## Workspace ##########################################################
# structure used for storing data for IPOPT solver
mutable struct IPOPTworkspace <: SubsolverWorkspace
    # original problem
    problem::Problem
	# user settings
	settings::IPOPTsettings
	# outdated flag
	outdated::Bool
	# sparsity info
	numVars::Int
	numCnss::Int
	numValsCnsJcb::Int
	numValsLgrHes::Int
	# problem functions
	evalCnsVal::Function
	evalCnsJcb::Function
	evalObjVal::Function
	evalObjGrd::Function
	evalLgrHes::Function
end


## Setup & Update ##########################################################
# this function creates an IPOPT.Model representing the given CvxQproblem
function setup(problem::Problem,settings::IPOPTsettings;withPrecompilation::Bool=false,nodeType::Type=AbstractNode)::IPOPTworkspace

	# check the problem
    @assert problem.objFun isa NullObjective || problem.objFun isa LinearObjective ||
			problem.objFun isa QuadraticObjective || problem.objFun isa ConvexObjective
    @assert problem.cnsSet isa NullConstraintSet || problem.cnsSet isa LinearConstraintSet ||
			problem.cnsSet isa ConvexConstraintSet

	# ensure type consistency
	varSet = problem.varSet
	if problem.cnsSet isa ConvexConstraintSet
		cnsSet = problem.cnsSet
	else
		cnsSet = ConvexConstraintSet{SpMatrix{Float},SpMatrix{Float}}(problem.cnsSet)
	end
	if problem.objFun isa ConvexObjective
		objFun = problem.objFun
	else
		objFun = ConvexObjective{SpMatrix{Float},Vector{Float}}(problem.objFun)
	end

	# collect info
	(numVars,numCnss) = (get_size(varSet),get_size(cnsSet))

	# build constraints jacobian sparsity info
	(rows_,cols_) = get_jacobianSparsity(cnsSet)
	(rowsCnsJcb,colsCnsJcb) = (@. Int32(rows_),@. Int32(cols_))
	numValsCnsJcb = length(rowsCnsJcb)

	# build the lagrangian hessian spasity info
	lgrHesSparsity = deepcopy(objFun.hesSparsity)
	for k in 1:numCnss
		(rows_,cols_,~) = findnz(cnsSet.hesSparsity[k])
		for i in 1:nnz(cnsSet.hesSparsity[k])
			lgrHesSparsity[rows_[i],cols_[i]] = true
		end
	end
	(rows_,cols_,~) = findnz(lgrHesSparsity)
	(rowsLgrHes,colsLgrHes) = (@. Int32(rows_),@. Int32(cols_))
	numValsLgrHes = length(rowsLgrHes)


	## build problem functions
	# objective and constraints values
	evalObjVal = objFun.evalVal
	evalCnsVal_ = cnsSet.evalVal
	function evalCnsVal(x::Vector{Float},cnsVals::VirtualVector{Float})::Nothing
		vals_ = evalCnsVal_(x)
		for k in 1:numCnss
			cnsVals[k] = vals_[k]
		end
		return
	end
	# objective gradient and constraints jacobian
	evalObjGrd_ = objFun.evalGrd
	function evalObjGrd(x::Vector{Float},objGrad::Vector{Float})::Nothing
		grd_ = evalObjGrd_(x)
		for k in 1:numVars
			objGrad[k] = grd_[k]
		end
		return
	end
	evalCnsJcb_ = cnsSet.evalJcb
	function evalCnsJcb(x::Vector{Float},mode::Symbol,rows::VirtualVector{Int32},cols::VirtualVector{Int32},vals::VirtualVector{Float})::Nothing
		if mode == :Structure
			for k in 1:numValsCnsJcb
				rows[k] = rowsCnsJcb[k]
				cols[k] = colsCnsJcb[k]
			end
		elseif mode == :Values
			jacobian_ = evalCnsJcb_(x)
			for k in 1:numValsCnsJcb
				vals[k] = jacobian_[rowsCnsJcb[k],colsCnsJcb[k]]
			end
		end
		return
	end
	# lagrangian hessian
	evalObjHes_ = objFun.evalHes
	evalCnsHes_ = cnsSet.evalHes
	function evalLgrHes(x::Vector{Float},mode::Symbol,rows::Vector{Int32},cols::Vector{Int32},objFactor::Float,lambda::VirtualVector{Float},vals::Vector{Float})::Nothing
		if mode == :Structure
			for k in 1:numValsLgrHes
				rows[k] = rowsLgrHes[k]
				cols[k] = colsLgrHes[k]
			end
		elseif mode == :Values
			hessian_ = objFactor*evalObjHes_(x)
			if numCnss>0
				hessian_ += sum(lambda.*evalCnsHes_(x))
			end
			for k in 1:numValsLgrHes
				vals[k] = hessian_[rowsLgrHes[k],colsLgrHes[k]]
			end
		end
		return
	end



	# create the workspace
	workspace = IPOPTworkspace(problem,settings,false,
							   numVars,numCnss,
							   numValsCnsJcb,numValsLgrHes,
							   evalCnsVal,evalCnsJcb,
							   evalObjVal,evalObjGrd,
							   evalLgrHes)

	# precompile the main functions to be used according to the workspace created
	if withPrecompilation
		precompile(make_outdated!,(typeof(workspace),))
		precompile(update!,(typeof(workspace),))
		precompile(solve!,(nodeType,typeof(workspace)))
	end


	return workspace


end

# it marks the workspace as outdated
function make_outdated!(workspace::IPOPTworkspace)::Nothing
    workspace.outdated = true
    return
end

#
function update!(workspace::IPOPTworkspace)::Nothing

	# create a new workspace (the interface lacks of modifying functions)
	workspace_ = setup(workspace.problem,workspace.settings)

	# copy the data contained in the new workspace into the old
	workspace.numVars = workspace_.numVars
	workspace.numCnss = workspace_.numCnss
	workspace.numValsCnsJcb = workspace_.numValsCnsJcb
	workspace.numValsLgrHes = workspace_.numValsLgrHes
	workspace.evalCnsVal = workspace_.evalCnsVal
	workspace.evalCnsJcb = workspace_.evalCnsJcb
	workspace.evalObjVal = workspace_.evalObjVal
	workspace.evalObjGrd = workspace_.evalObjGrd
	workspace.evalLgrHes = workspace_.evalLgrHes

    # mark the workspace as up to date
    workspace.outdated = false

    return
end

## Solve ##########################################################
function solve!(node::AbstractNode,workspace::IPOPTworkspace;objUpperLimit::Float=Inf,withInfeasibilityMinimization::Bool=false)::Tuple{Int8,Float}

	# check if the problem is outdated and update it
	if workspace.outdated
		update!(workspace)
	end

	# collect info on the problem
	numVars = workspace.numVars
	numCnss = workspace.numCnss
	withCuts = nnz(SpMatrix{Float}(node.cutSet.A)) > 0

	# include the possible local cuts in the formulation
	if withCuts
		cutSet = ConvexConstraintSet{SpMatrix{Float},SpMatrix{Float}}(node.cutSet)
		numCuts = get_size(cutSet)

		# build new constraint evaluation function
		oldEvalCnsVal = workspace.evalCnsVal
		function evalCnsCutVal(x::Vector{Float},val::Vector{Float})::Nothing
			@assert length(val) >= numCuts
			oldEvalCnsVal(x,@view val[1:numCnss])
			val[numCnss+1:numCnss+numCuts] .= cutSet.evalVal(x)
			return
		end

		# build new jacobian function (with sparsity info)
		offset_ = Int32(numCnss)
		(rows_,cols_) = get_jacobianSparsity(cutSet)
		(rowsCutJcb,colsCutJcb) = (@. Int32(rows_),@. Int32(cols_))
		numValsCutJcb = length(rowsCutJcb)
		function evalCutJcb(x::Vector{Float},mode::Symbol,rows::VirtualVector{Int32},cols::VirtualVector{Int32},vals::VirtualVector{Float})::Nothing
			if mode == :Structure
				@assert length(rows) == length(cols) >= numValsCutJcb
				@. rows = rowsCutJcb + offset_
				@. cols = colsCutJcb
			elseif mode == :Values
				@assert length(vals) >= numValsCutJcb
				jacobian_ = cutSet.evalJcb(x)
				@. vals = getindex([jacobian_],rowsCutJcb,colsCutJcb)
			end
			return
		end

		# build new constraints jacobian function
		numValsCnsJcb = workspace.numValsCnsJcb
		oldEvalCnsJcb = workspace.evalCnsJcb

		function evalCnsCutJcb(x::Vector{Float},mode::Symbol,rows::Vector{Int32},cols::Vector{Int32},vals::Vector{Float})::Nothing
			oldEvalCnsJcb(x,mode,view(rows,1:numValsCnsJcb),view(cols,1:numValsCnsJcb),view(vals,1:numValsCnsJcb))
			evalCutJcb(x,mode,view(rows,numValsCnsJcb+1:numValsCnsJcb+numValsCutJcb),view(cols,numValsCnsJcb+1:numValsCnsJcb+numValsCutJcb),view(vals,numValsCnsJcb+1:numValsCnsJcb+numValsCutJcb))
			return
		end

		# wrap the hessian function
		oldEvalLgrHes = workspace.evalLgrHes
		function evalLgrHes(x::Vector{Float},mode::Symbol,rows::Vector{Int32},cols::Vector{Int32},objFactor::Float,lambda::Vector{Float},vals::Vector{Float})::Nothing
			oldEvalLgrHes(x,mode,rows,cols,objFactor,view(lambda,1:numCnss),vals)
			return
		end

		# build an ipopt problem model
		model = createProblem(numVars, node.varLoBs, node.varUpBs,
							  numCnss+numCuts, vcat(node.cnsLoBs,cutSet.loBs),
							  				   vcat(node.cnsUpBs,cutSet.upBs),
							  workspace.numValsCnsJcb+numValsCutJcb, workspace.numValsLgrHes,
							  workspace.evalObjVal, evalCnsCutVal,
							  workspace.evalObjGrd, evalCnsCutJcb,
							  evalLgrHes)


	else

		# build an ipopt problem model
		model = createProblem(numVars,node.varLoBs,node.varUpBs,
							  numCnss,node.cnsLoBs,node.cnsUpBs,
							  workspace.numValsCnsJcb,workspace.numValsLgrHes,
							  workspace.evalObjVal,workspace.evalCnsVal,
							  workspace.evalObjGrd,workspace.evalCnsJcb,
							  workspace.evalLgrHes)


	end

	# default settings for OpenBB, taking care of the fact that some defaults in OpenBB are different from Ipopt defaults
	defaultSettings = IPOPTsettings(print_level=5,expect_infeasible_problem="no",mu_strategy="monotone",mu_oracle="quality-function",bound_relax_factor=1e-8)
	# set non-default settings
	for field in fieldnames(IPOPTsettings)
		if getfield(workspace.settings,field) != getfield(defaultSettings,field)
			addOption(model,String(field),getfield(workspace.settings,field))
		end
	end

	# finally solve the problem
	runtime = @elapsed ipoptStatus = solveProblem(model)
	# collect results
	if ipoptStatus in [0,1]
		status = 0 # "solved"
		# store the solution in the node:
		@. node.primal = model.x
		@. node.bndDual = model.mult_x_U - model.mult_x_L
		@. node.cnsDual = model.mult_g[1:numCnss]
		if withCuts
			@. node.cutDual = model.mult_g[numCnss+1:numCnss+numCuts]
		end
		node.objUpB = model.obj_val
		node.objLoB = node.objUpB # assume perfect solution

	elseif ipoptStatus in [2]
		status = 1 # "infeasible"

		# store the solution in the node:
		@. node.primal = model.x
		@. node.bndDual = model.mult_x_U - model.mult_x_L
		@. node.cnsDual = model.mult_g[1:numCnss]
		if withCuts
			@. node.cutDual = model.mult_g[numCnss+1:numCnss+numCuts]
		end
		node.objUpB = Inf
		node.objLoB = Inf

	elseif ipoptStatus in [3,4,6,-1,-2,-4]
		status = 2 # unreliable
		# store the solution in the node:
		@. node.primal = model.x
		@. node.bndDual = model.mult_x_U - model.mult_x_L
		@. node.cnsDual = model.mult_g[1:numCnss]
		if withCuts
			@. node.cutDual = model.mult_g[numCnss+1:numCnss+numCuts]
		end
		# compute objective lower and upper bounds
		oldObjLoB = node.objLoB
		update_objBounds!(node,workspace.problem,
						  workspace.settings.acceptable_constr_viol_tol,
					      workspace.settings.acceptable_dual_inf_tol)
		node.objLoB = max(node.objLoB,oldObjLoB) # it is possible that the parent node solution was more successful
		@warn "Ipopt, unreliable solve: status = "*string(Ipopt.ApplicationReturnStatus[ipoptStatus],", code = ",ipoptStatus)
	else
		error("Ipopt "*string(Ipopt.ApplicationReturnStatus[ipoptStatus]," code = ",ipoptStatus))
	end

    return (status,runtime)
end
