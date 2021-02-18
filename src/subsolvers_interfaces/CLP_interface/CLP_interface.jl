# @Author: Massimo De Mauri <massimo>
# @Date:   2020-01-08T14:39:40+01:00
# @Email:  massimo.demauri@gmail.com
# @Filename: CLP_interface.jl
# @Last modified by:   massimo
# @Last modified time: 2021-02-11T17:43:34+01:00
# @License: LGPL-3.0
# @Copyright: {{copyright}}

try CCLP
catch err include("./CLP_Cwrapper.jl") end



mutable struct CLPsettings <: SubsolverSettings
	algorithm::Int
	primalTolerance::Float64
	dualTolerance::Float64
	objectiveOffset::Float64
	dualObjectiveLimit::Float64
	maximumIterations::Int
	maximumSeconds::Float64
	logLevel::Int
	scalingMode::Int
	perturbationMode::Int
end

function CLPsettings(;  algorithm::Int=1,
						primalTolerance::Float64=1e-6,
                        dualTolerance::Float64=1e-6,
						objectiveOffset::Float64=0.0,
                        dualObjectiveLimit::Float64=1.7976931348623157e308,
                        maximumIterations::Int=100000,
                        maximumSeconds::Float64=-1.0,
                        logLevel::Int=0,
                        scalingMode::Int=0,
                        perturbationMode::Int=100)::CLPsettings

   return CLPsettings(algorithm,primalTolerance,dualTolerance,objectiveOffset,
   					  dualObjectiveLimit,maximumIterations,maximumSeconds,
					  logLevel,scalingMode,perturbationMode)::CLPsettings
end

## Utility Functions ##########################################################

# wrapper to get verbosity of each solver in the same way
function get_verbosity(settings::CLPsettings)::Bool
	return settings.logLevel > 0
end

# wrapper to get primal tolerance of each solver in the same way
function get_primalTolerance(settings::CLPsettings)::Float
	return settings.primalTolerance
end

# wrapper to get dual tolerance of each solver in the same way
function get_dualTolerance(settings::CLPsettings)::Float
	return settings.dualTolerance
end

# wrapper to get iterations limit of each solver in the same way
function get_iterationsLimit(settings::CLPsettings)::Int
	return settings.maximumIterations
end

# wrapper to get time limit of each solver in the same way
function get_timeLimit(settings::CLPsettings)::Float
	return settings.maximumSeconds
end

# wrapper to set verbosity of each solver in a consistent way
function set_verbosity!(settings::CLPsettings,verbose::Bool)::Nothing
	if verbose
		settings.logLevel = 1
	end
	return
end

# wrapper to set primal tolerance of each solve in a consistent way
function set_primalTolerance!(settings::CLPsettings,tolerance::Float)::Nothing
	settings.primalTolerance = tolerance
	return
end

# wrapper to set dual tolerance of each solve in a consistent way
function set_dualTolerance!(settings::CLPsettings,tolerance::Float)::Nothing
	settings.dualTolerance = tolerance
	return
end

# wrapper to set the iteration limit of each solver in a consistent way
function set_iterationsLimit!(settings::CLPsettings,limit::Int)::Nothing
	# set the appropriate setting to the value
	settings.maximumIterations=limit
	return
end

# wrapper to set the time limit of each solver in a consistent way
function set_timeLimit!(settings::CLPsettings,limit::Float)::Nothing
	# set the appropriate setting to the value
	settings.maximumSeconds=limit
	return
end


mutable struct CLPworkspace <: SubsolverWorkspace
   # problem
   problem::Problem
   # memory
   prbData::Tuple{SpMatrix{Float},Vector{Float},SpMatrix{Float}}
   settings::CLPsettings
   # outdated flag
   outdated::Bool
end


## Setup & Update ##########################################################
# this function creates an Clp model representing the given linear problem
function setup(problem::Problem,settings::CLPsettings;withPrecompilation::Bool=false,nodeType::Type=AbstractNode)::CLPworkspace

   # check the problem
   @assert problem.objFun isa NullObjective || problem.objFun isa LinearObjective || problem.objFun isa QuadraticObjective
   @assert problem.cnsSet isa NullConstraintSet || problem.cnsSet isa LinearConstraintSet


	# ensure type consistency
	varSet = problem.varSet
	if problem.objFun isa QuadraticObjective
		objFun = QuadraticObjective{SpMatrix{Float},Vector{Float}}(problem.objFun)
	elseif problem.objFun isa LinearObjective
		objFun = LinearObjective{Vector{Float}}(problem.objFun)
	else
		objFun = NullObjective()
	end
	cnsSet = LinearConstraintSet{SpMatrix{Float}}(problem.cnsSet)


	# collect info
	numVars = get_size(problem.cnsSet)
	if objFun isa QuadraticObjective
		prbData = (cnsSet.A,objFun.L,objFun.Q)
	else
		prbData = (cnsSet.A,objFun.L,spzeros(numVars,numVars))
	end

	# create the Clp workspace
    workspace = CLPworkspace(problem,prbData,settings,false)

	# precompile the main functions to be used according to the workspace created
	if withPrecompilation
		precompile(make_outdated!,(typeof(workspace),))
		precompile(update!,(typeof(workspace),))
		precompile(solve!,(nodeType,typeof(workspace)))
	end

    # create the Clp workspace
    return workspace
end


# it marks the workspace as outdated
function make_outdated!(workspace::CLPworkspace)::Nothing
    workspace.outdated = true
    return
end


#
function update!(workspace::CLPworkspace)::Nothing
	# setup again
	workspace.prbData = setup(workspace.problem,workspace.settings,withPrecompilation=false).prbData
    # mark the workspace as up to date
    workspace.outdated = false
    return
end


## Solve ##########################################################
function solve!(node::BBnode,workspace::CLPworkspace;objUpperLimit::Float=Inf)::Tuple{Int8,Float64}

	# update the workspace if necessary
	if workspace.outdated
		update!(workspace)
	end

   # collect info on the problem
	numVars = get_size(workspace.problem.varSet)
	numCnss = get_size(workspace.problem.cnsSet)
	withCuts = nnz(sparse(node.cutSet.A)) > 0

	_model = CCLP.ClpModel()
	update_settings!(_model,workspace.settings)
	if !withCuts
		CCLP.load_problem(_model,workspace.prbData[2],node.varLoBs,node.varUpBs,
								 workspace.prbData[1],node.cnsLoBs,node.cnsUpBs)
		if nnz(workspace.prbData[3]) > 0
			CCLP.load_quadratic_objective(_model,workspace.prbData[3])
		end
	else
		CCLP.load_problem(_model,workspace.prbData[2],node.varLoBs,node.varUpBs,
						  vcat(workspace.prbData[1],node.cutSet.A),vcat(node.cnsLoBs,node.cutSet.loBs),vcat(node.cnsUpBs,node.cutSet.upBs))
		if nnz(workspace.prbData[3]) > 0
			CCLP.load_quadratic_objective(_model,workspace.prbData[3])
		end
	end

	# set dual limit
	# CCLP.set_dual_objective_limit(_model,objUpperLimit - workspace.problem.objFun.c)

	# set hotstart info
	#TODO

	# solve problem
	if workspace.settings.algorithm == 0 # dual
		CLPtime = @elapsed CLPstatus = CCLP.dual(_model,1)
	elseif workspace.settings.algorithm == 1 # primal
		CLPtime = @elapsed CLPstatus = CCLP.primal(_model,1)
	elseif workspace.settings.algorithm == 3 # barrier
		CLPtime = @elapsed CLPstatus = CCLP.initial_barrier_solve(_model)
	else
		@error "CLP : algorithm unavailable"
	end

	# copy results in the node
	primal = CCLP.primal_column_solution(_model)[1:numVars]
	node.primal .= primal
	node.bndDual .= -CCLP.get_reduced_cost(_model)[1:numVars]
	node.cnsDual .= -CCLP.get_row_price(_model)[1:numCnss]
	if withCuts
		node.cutDual .= -CCLP.get_row_price(_model)[numCnss+1:numCnss+numCuts]
	end

	if CLPstatus == 0
		status = 0 # solved
		update_objBounds!(node,workspace.problem,10*workspace.settings.primalTolerance,10*workspace.settings.dualTolerance)
		node.reliable = true
	elseif CLPstatus in [1,2]
       status = 1 # infeasible or unbounded
	   node.objUpB = Inf
	   node.objLoB = Inf
	   node.reliable = true
   	elseif CLPstatus == 3
		update_objBounds!(node,workspace.problem,10*workspace.settings.primalTolerance,10*workspace.settings.dualTolerance)
		if CCLP.is_iteration_limit_reached(_model) ||
		   CLPtime >= workspace.settings.maximumSeconds ||
		   CCLP.is_dual_objective_limit_reached(_model)
			status = 1 # solved
			reliable = true
		else
			status = 2 # unreliable
			@warn "Inaccuracy in node sol, status code: "*string(CLPstatus)
			node.reliable = false
		end
   	elseif CLPstatus == 4
       status = 3 # "error"
        @error "Subsover error, status code: "*string(CLPstatus)
	end

    return (status,CLPtime)
end


# updates the CLP settings
function update_settings!(model::CCLP.ClpModel,settings::CLPsettings)::Nothing
	CCLP.set_primal_tolerance(model,settings.primalTolerance)
	CCLP.set_dual_tolerance(model,settings.dualTolerance)
	CCLP.set_objective_offset(model,settings.objectiveOffset)
	CCLP.set_dual_objective_limit(model,settings.dualObjectiveLimit)
	CCLP.set_maximum_iterations(model,settings.maximumIterations)
	CCLP.set_maximum_seconds(model,settings.maximumSeconds)
	CCLP.set_log_level(model,settings.logLevel)
	CCLP.scaling(model,settings.scalingMode)
	CCLP.set_perturbation(model,settings.perturbationMode)
	return
end
