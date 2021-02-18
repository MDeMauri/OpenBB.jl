# @Author: Massimo De Mauri <massimo>
# @Date:   2021-01-27T18:47:46+01:00
# @Email:  massimo.demauri@protonmail.com
# @Filename: JuMP.jl
# @Last modified by:   massimo
# @Last modified time: 2021-02-09T20:02:28+01:00
# @License: LGPL-3.0
# @Copyright: {{copyright}}
using JuMP
using CPLEX

const JuMP_Inf = 1e20

mutable struct JuMP_CPLEXsettings <:SubsolverSettings
	Simplex_Display::Bool
	ParamDisplay::Bool
    Preprocessing_Presolve::Bool
    LPMethod::Int
	QPMethod::Int
	SolutionType::Int
	#
	PrimalTolerance::Float
	DualTolerance::Float
    IterationsLimit::Int
	# algorithm modifiers
	Threads::Int
	TimeLimit::Float
end

function JuMP_CPLEXsettings(;
						Simplex_Display::Bool=false,
						ParamDisplay::Bool=false,
						Preprocessing_Presolve::Bool=true,
                        LPMethod::Int=2,
						QPMethod::Int=2,
						SolutionType::Int=2,
						PrimalTolerance::Float=1e-4,
						DualTolerance::Float=1e-4,
						IterationsLimit::Int=9223372036800000000,
						Threads::Int=1,
						TimeLimit::Float=1e20)::JuMP_CPLEXsettings
    return JuMP_CPLEXsettings(
						 Simplex_Display,
						 ParamDisplay,
						 Preprocessing_Presolve,
						 LPMethod,
						 QPMethod,
						 SolutionType,
						 PrimalTolerance,
						 DualTolerance,
						 IterationsLimit,
						 Threads,
						 TimeLimit)
end





## Utility Functions ##########################################################


# wrapper to get verbosity of each solver in the same way
function get_verbosity(settings::JuMP_CPLEXsettings)::Bool
	return false #TODO
end

# wrapper to get primal tolerance of each solver in the same way
function get_primalTolerance(settings::JuMP_CPLEXsettings)::Float
	return settings.PrimalTolerance
end

# wrapper to get dual tolerance of each solver in the same way
function get_dualTolerance(settings::JuMP_CPLEXsettings)::Float
	return settings.DualTolerance
end

# wrapper to get iterations limit of each solver in the same way
function get_iterationsLimit(settings::JuMP_CPLEXsettings)::Int
	return settings.IterationsLimit
end

# wrapper to get time limit of each solver in the same way
function get_timeLimit(settings::JuMP_CPLEXsettings)::Float
	return ettings.TimeLimit
end

# wrapper to set verbosity of each solver in a consistent way
function set_verbosity!(settings::JuMP_CPLEXsettings,verbose::Bool)::Nothing
	if verbose
		#TODO
	end
	return
end

# wrapper to set primal tolerance of each solve in a consistent way
function set_primalTolerance!(settings::JuMP_CPLEXsettings,tolerance::Float)::Nothing
	settings.PrimalTolerance = tolerance
	return
end

# wrapper to set dual tolerance of each solve in a consistent way
function set_dualTolerance!(settings::JuMP_CPLEXsettings,tolerance::Float)::Nothing
	settings.DualTolerance = tolerance
	return
end

# wrapper to set the iteration limit of each solver in a consistent way
function set_iterationsLimit!(settings::JuMP_CPLEXsettings,limit::Int)::Nothing
	# set the appropriate setting to the value
	settings.IterationsLimit=limit
	return
end

# wrapper to set the time limit of each solver in a consistent way
function set_timeLimit!(settings::JuMP_CPLEXsettings,limit::Float)::Nothing
	# set the appropriate setting to the value
	settings.TimeLimit=limit
	return
end



## Workspace ##########################################################
# structure used for storing data for Clp solver
mutable struct JuMP_CPLEXworkspace <: SubsolverWorkspace
   # problem
   problem::Problem
   # memory
   model::JuMP.Model
   # settings
   settings::JuMP_CPLEXsettings
   # outdated flag
   outdated::Bool
end


## Setup & Update ##########################################################
# this function creates an Clp model representing the given linear problem
function setup(problem::Problem,settings::JuMP_CPLEXsettings;withPrecompilation::Bool=false,nodeType::Type=AbstractNode)::JuMP_CPLEXworkspace

   # check the problem
   @assert problem.objFun isa NullObjective || problem.objFun isa LinearObjective || problem.objFun isa QuadraticObjective
   @assert problem.cnsSet isa NullConstraintSet || problem.cnsSet isa LinearConstraintSet

	# create a CPLEX model
	model = Model(CPLEX.Optimizer); update_settings!(model,settings)

	# collect info
	numVars = get_size(problem.varSet)
	numCnss = get_size(problem.cnsSet)

	# define problem
	@variable(model,problem.varSet.loBs[j] <= vars[j=1:numVars] <= problem.varSet.upBs[j])
	@variable(model,problem.cnsSet.loBs[j] <= slacks[j=1:numCnss] <= problem.cnsSet.upBs[j])
	@constraint(model,base,problem.cnsSet.A*vars - slacks .== 0)
	if problem.objFun isa QuadraticObjective
		@objective(model,Min,.5*vars'*problem.objFun.Q*vars + problem.objFun.L'*vars)
	elseif problem.objFun isa LinearObjective
		@objective(model,Min,problem.objFun.L'*vars)
	end

	# create the Clp workspace
	workspace = JuMP_CPLEXworkspace(problem,model,settings,false)

	# precompile the main functions to be used according to the workspace created
	if withPrecompilation
		@assert precompile(make_outdated!,(typeof(workspace),))
		@assert precompile(update!,(typeof(workspace),))
		@assert precompile(solve!,(nodeType,typeof(workspace)))
	end

    return workspace
end


# it marks the workspace as outdated
function make_outdated!(workspace::JuMP_CPLEXworkspace)::Nothing
    workspace.outdated = true
    return
end


function update!(workspace::JuMP_CPLEXworkspace)::Nothing

	# check the problem
   @assert workspace.problem.objFun isa NullObjective || workspace.problem.objFun isa LinearObjective || workspace.problem.objFun isa QuadraticObjective
   @assert workspace.problem.cnsSet isa NullConstraintSet || workspace.problem.cnsSet isa LinearConstraintSet

	# setup again
	workspace.model = setup(workspace.problem,workspace.settings).model

    # mark the workspace as up to date
    workspace.outdated = false
    return
end


## Solve ##########################################################
function solve!(node::AbstractNode,workspace::JuMP_CPLEXworkspace;objUpperLimit::Float=Inf)::Tuple{Int8,Float}

	# keep timing
	start_time = time()

	# update the workspace if necessary
	if workspace.outdated
		update!(workspace)
	end

   # collect info on the problem
	numVars = get_size(workspace.problem.varSet)
	numCnss = get_size(workspace.problem.cnsSet)
	numCuts = get_size(node.cutSet)
	withCuts = nnz(node.cutSet.A) > 0

	# adapt the problem definition to the node
	set_lower_bound.(object_dictionary(workspace.model)[:vars],node.varLoBs)
	set_upper_bound.(object_dictionary(workspace.model)[:vars],node.varUpBs)
	set_lower_bound.(object_dictionary(workspace.model)[:slacks],node.cnsLoBs)
	set_upper_bound.(object_dictionary(workspace.model)[:slacks],node.cnsUpBs)

	if withCuts # there are some local cuts
		@variable(workspace.model,node.cutSet.loBs <= cutSlacks[j=1:numCuts] <= node.cutSet.upBs)
		@constraint(workspace.model,cuts,node.cutSet.A .== cutSlacks)
	end

	# set the objective limit
	# set_optimizer_attribute(workspace.model,"CPXPARAM_Simplex_Limits_UpperObj",objUpperLimit - workspace.problem.objFun.c)

	# set hotstart info
	#TODO: check if this is possible

	# solve
	optimize!(workspace.model)
	jumpStatus = Symbol(termination_status(workspace.model))

	# collect solution
	try @. node.primal  = value(object_dictionary(workspace.model)[:vars]) catch e end
	try @. node.bndDual = -dual(UpperBoundRef(object_dictionary(workspace.model)[:vars])) +
						  -dual(LowerBoundRef(object_dictionary(workspace.model)[:vars])) catch e end
	try @. node.cnsDual = -dual(object_dictionary(workspace.model)[:base]) catch e end
	if withCuts
		try @. node.cutDual = -dual(object_dictionary(workspace.model)[:cuts]) catch e end
	end


	if jumpStatus in [:OPTIMAL,:LOCALLY_SOLVED]
		status = 0
		update_objBounds!(node,workspace.problem,10*workspace.settings.PrimalTolerance,10*workspace.settings.DualTolerance)
	elseif jumpStatus in [:OBJECTIVE_LIMIT]
		status = 0
		(node.objLoB,node.objUpB) = (objUpperLimit,Inf)
	elseif jumpStatus in [:INFEASIBLE,:DUAL_INFEASIBLE,:LOCALLY_INFEASIBLE,:INFEASIBLE_OR_UNBOUNDED]
		status = 1
		(node.objLoB,node.objUpB) = (Inf,Inf)
	elseif jumpStatus in [:ALMOST_OPTIMAL,:ALMOST_INFEASIBLE,:ALMOST_DUAL_INFEASIBLE,:ALMOST_LOCALLY_SOLVED,:SLOW_PROGRESS,
						  :ITERATION_LIMIT,:TIME_LIMIT,:NODE_LIMIT,:SOLUTION_LIMIT,:MEMORY_LIMIT,:NORM_LIMIT]
		status = 2
		update_objBounds!(node,workspace.problem,10*workspace.settings.PrimalTolerance,10*workspace.settings.DualTolerance)
		@warn "Inaccuracy in node sol, status : "*string(jumpStatus)
	elseif jumpStatus in [:OPTIMIZE_NOT_CALLED,:NUMERICAL_ERROR,:INVALID_MODEL,:INVALID_OPTION,:INTERRUPTED,:OTHER_ERROR,:OTHER_LIMIT]
		status = 3
		error("Subsover error, status: "*string(jumpStatus))
	end

	# clean up
	if withCuts
		delete(workspace.model,object_dictionary(workspace.model)[:cutSlacks])
		delete(workspace.model,object_dictionary(workspace.model)[:cuts])
	end

    return (status,time()-start_time)
end


## helper functions
function update_settings!(model::JuMP.Model,settings::JuMP_CPLEXsettings)::Nothing
	# communicate settings to JuMP
	settings_dict = Dict{String,Any}("CPXPARAM_Simplex_Tolerances_Feasibility"=>0.5*settings.PrimalTolerance,
									 "CPXPARAM_Simplex_Tolerances_Optimality"=>0.5*settings.DualTolerance,
									 "CPXPARAM_Simplex_Limits_Iterations"=>settings.IterationsLimit)


	for field in setdiff(fieldnames(JuMP_CPLEXsettings),[:PrimalTolerance,:DualTolerance,:IterationsLimit])
		settings_dict["CPXPARAM_"*String(field)] = getfield(settings,field)
	end
	set_optimizer_attributes(model,settings_dict...)
	return
end
