# @Author: Massimo De Mauri <massimo>
# @Date:   2021-01-27T18:47:46+01:00
# @Email:  massimo.demauri@protonmail.com
# @Filename: JuMP.jl
# @Last modified by:   massimo
# @Last modified time: 2021-02-09T20:02:36+01:00
# @License: LGPL-3.0
# @Copyright: {{copyright}}
using JuMP
using Clp

const JuMP_Inf = 1e20

mutable struct JuMP_CLPsettings <: SubsolverSettings
	SolveType::Int
	PrimalTolerance::Float
	DualTolerance::Float
	DualObjectiveLimit::Float
	MaximumIterations::Int
	MaximumSeconds::Float
	LogLevel::Int
	PresolveType::Int
	InfeasibleReturn::Int
	Scaling::Int
	Perturbation::Int
end

function JuMP_CLPsettings(; SolveType::Int=0,
					        PrimalTolerance::Float=1e-4, # actual solver default: 1e-7
	                        DualTolerance::Float=1e-6, # actual solver default: 1e-7
	                        DualObjectiveLimit::Float=1.7976931348623157e308,
	                        MaximumIterations::Int=100000,
	                        MaximumSeconds::Float=-1.0,
	                        LogLevel::Int=0, # actual solver default: 1
							PresolveType::Int=0,
							InfeasibleReturn::Int=1, # actual solver default: 0
	                        Scaling::Int=3,
	                        Perturbation::Int=100)::JuMP_CLPsettings

   return JuMP_CLPsettings(SolveType,PrimalTolerance,DualTolerance,DualObjectiveLimit,
   					  	   MaximumIterations,MaximumSeconds,LogLevel,PresolveType,
					  	   InfeasibleReturn,Scaling,Perturbation)::JuMP_CLPsettings
end




## Utility Functions ##########################################################

# wrapper to get verbosity of each solver in the same way
function get_verbosity(settings::JuMP_CLPsettings)::Bool
	return settings.LogLevel > 0
end

# wrapper to get primal tolerance of each solver in the same way
function get_primalTolerance(settings::JuMP_CLPsettings)::Float
	return settings.PrimalTolerance
end

# wrapper to get dual tolerance of each solver in the same way
function get_dualTolerance(settings::JuMP_CLPsettings)::Float
	return settings.DualTolerance
end

# wrapper to get iterations limit of each solver in the same way
function get_iterationsLimit(settings::JuMP_CLPsettings)::Int
	return settings.MaximumIterations
end

# wrapper to get time limit of each solver in the same way
function get_timeLimit(settings::JuMP_CLPsettings)::Float
	return settings.MaximumSeconds
end

# wrapper to set verbosity of each solver in a consistent way
function set_verbosity!(settings::JuMP_CLPsettings,verbose::Bool)::Nothing
	if verbose
		settings.LogLevel = 1
	end
	return
end

# wrapper to set primal tolerance of each solve in a consistent way
function set_primalTolerance!(settings::JuMP_CLPsettings,tolerance::Float)::Nothing
	settings.PrimalTolerance = tolerance
	return
end

# wrapper to set dual tolerance of each solve in a consistent way
function set_dualTolerance!(settings::JuMP_CLPsettings,tolerance::Float)::Nothing
	settings.DualTolerance = tolerance
	return
end

# wrapper to set the iteration limit of each solver in a consistent way
function set_iterationsLimit!(settings::JuMP_CLPsettings,limit::Int)::Nothing
	# set the appropriate setting to the value
	settings.MaximumIterations=limit
	return
end

# wrapper to set the time limit of each solver in a consistent way
function set_timeLimit!(settings::JuMP_CLPsettings,limit::Float)::Nothing
	# set the appropriate setting to the value
	settings.MaximumSeconds=limit
	return
end


## Workspace ##########################################################
# structure used for storing data for Clp solver
mutable struct JuMP_CLPworkspace <: SubsolverWorkspace
   # problem
   problem::Problem
   # memory
   model::JuMP.Model
   # settings
   settings::JuMP_CLPsettings
   # outdated flag
   outdated::Bool
end


## Setup & Update ##########################################################
# this function creates an Clp model representing the given linear problem
function setup(problem::Problem,settings::JuMP_CLPsettings;withPrecompilation::Bool=false,nodeType::Type=AbstractNode)::JuMP_CLPworkspace

   # check the problem
   @assert problem.objFun isa NullObjective || problem.objFun isa LinearObjective
   @assert problem.cnsSet isa NullConstraintSet || problem.cnsSet isa LinearConstraintSet

	# create a Clp model
	model = Model(Clp.Optimizer);update_settings!(model,settings)

	# collect info
	numVars = get_size(problem.varSet)
	numCnss = get_size(problem.cnsSet)

	# define problem
	@variable(model,problem.varSet.loBs[j] <= vars[j=1:numVars] <= problem.varSet.upBs[j])
	@variable(model,problem.cnsSet.loBs[j] <= slacks[j=1:numCnss] <= problem.cnsSet.upBs[j])
	@constraint(model,base,problem.cnsSet.A*vars - slacks .== 0)
	if problem.objFun isa LinearObjective
		@objective(model,Min,problem.objFun.L'*vars)
	end

	# create the Clp workspace
	workspace = JuMP_CLPworkspace(problem,model,settings,false)

	# precompile the main functions to be used according to the workspace created
	if withPrecompilation
		@assert precompile(make_outdated!,(typeof(workspace),))
		@assert precompile(update!,(typeof(workspace),))
		@assert precompile(solve!,(nodeType,typeof(workspace)))
	end

    return workspace
end


# it marks the workspace as outdated
function make_outdated!(workspace::JuMP_CLPworkspace)::Nothing
    workspace.outdated = true
    return
end


function update!(workspace::JuMP_CLPworkspace)::Nothing

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
function solve!(node::AbstractNode,workspace::JuMP_CLPworkspace;objUpperLimit::Float=Inf)::Tuple{Int8,Float}

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
	# set_optimizer_attribute(workspace.model,"DualObjectiveLimit",objUpperLimit - workspace.problem.objFun.c)

	# set hotstart info
	#TODO: check if this is possible

	# solve
	optimize!(workspace.model)
	jumpStatus = Symbol(termination_status(workspace.model))

	# collect solution
	try @. node.primal  = value(object_dictionary(workspace.model)[:vars]) catch e end
	try @. node.bndDual = -dual(UpperBoundRef(object_dictionary(workspace.model)[:vars])) + -dual(LowerBoundRef(object_dictionary(workspace.model)[:vars])) catch e end
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


function update_settings!(model::JuMP.Model,settings::JuMP_CLPsettings)::Nothing

	# communicate settings to the solver
	settings_dict = Dict{String,Any}("PrimalTolerance"=>0.5*settings.PrimalTolerance,
									 "DualTolerance"=>0.5*settings.DualTolerance)
	for field in setdiff(fieldnames(JuMP_CLPsettings),[:PrimalTolerance,:DualTolerance])
		settings_dict[String(field)] = getfield(settings,field)
	end
	set_optimizer_attributes(model,settings_dict...)
	return
end
