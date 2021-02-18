# @Author: Massimo De Mauri <massimo>
# @Date:   2019-02-19T20:29:29+01:00
# @Email:  massimo.demauri@gmail.com
# @Filename: QPALM_interface.jl
# @Last modified by:   massimo
# @Last modified time: 2021-02-12T17:38:09+01:00
# @License: LGPL-3.0
# @Copyright: {{copyright}}

using QPALM
const maxQPALMiterations = 20000

## Settings ##########################################################

# structure used to hold the settings for QPALM
mutable struct QPALMsettings <: SubsolverSettings
	max_iter::Int64
	inner_max_iter::Int64
	eps_abs::Float
	eps_rel::Float
	eps_abs_in::Float
	eps_rel_in::Float
	rho::Float
	eps_prim_inf::Float
	eps_dual_inf::Float
	theta::Float
	delta::Float
	sigma_max::Float
	proximal::Bool
	gamma_init::Float
	gamma_upd::Float
	gamma_max::Float
	scaling::Int64
	nonconvex::Bool
	verbose::Bool
	print_iter::Int64
	warm_start::Bool
	reset_newton_iter::Int64
	enable_dual_termination::Bool
	dual_objective_limit::Float
	time_limit::Float
end

function QPALMsettings(;max_iter::Int64=maxQPALMiterations,
						inner_max_iter::Int64=100,
						eps_abs::Float=1.0e-6,
						eps_rel::Float=1.0e-6,
						eps_abs_in::Float=1.0,
						eps_rel_in::Float=1.0,
						rho::Float=0.1,
						eps_prim_inf::Float=1.0e-4,
						eps_dual_inf::Float=1.0e-6,
						theta::Float=0.25,
						delta::Float=100.0,
						sigma_max::Float=1.0e9,
						proximal::Bool=true,
						gamma_init::Float=10.0,
						gamma_upd::Float=10.0,
						gamma_max::Float=1.0e7,
						scaling::Int64=2,
						nonconvex::Bool=false,
						verbose::Bool=false,
						print_iter::Int64=1,
						warm_start::Bool=true,
						reset_newton_iter::Int64=100,
						enable_dual_termination::Bool=false,
						dual_objective_limit::Float=1.0e20,
						time_limit::Float=1.0e20)::QPALMsettings



    return QPALMsettings(max_iter,inner_max_iter,eps_abs,eps_rel,eps_abs_in,eps_rel_in,
						 rho,eps_prim_inf,eps_dual_inf,theta,delta,sigma_max,proximal,
						 gamma_init,gamma_upd,gamma_max,scaling,nonconvex,verbose,print_iter,
						 warm_start,reset_newton_iter,enable_dual_termination,dual_objective_limit,
						 time_limit)
end


## Utility Functions ##########################################################

# wrapper to get verbosity of each solver in the same way
function get_verbosity(settings::QPALMsettings)::Bool
	return settings.verbose
end

# wrapper to get primal tolerance of each solver in the same way
function get_primalTolerance(settings::QPALMsettings)::Float
	return settings.eps_prim_inf
end

# wrapper to get dual tolerance of each solver in the same way
function get_dualTolerance(settings::QPALMsettings)::Float
	return settings.eps_dual_inf
end

# wrapper to get iterations limit of each solver in the same way
function get_iterationsLimit(settings::QPALMsettings)::Int
	return settings.max_iter
end

# wrapper to get time limit of each solver in the same way
function get_timeLimit(settings::QPALMsettings)::Float
	return settings.time_limit
end


# wrapper to set verbosity of each solver in a consistent way
function set_verbosity!(settings::QPALMsettings,verbose::Bool)::Nothing
	settings.verbose = verbose
	return
end

# wrapper to set primal tolerance of each solve in a consistent way
function set_primalTolerance!(settings::QPALMsettings,tolerance::Float)::Nothing
	settings.eps_prim_inf= tolerance
	return
end

# wrapper to set dual tolerance of each solve in a consistent way
function set_dualTolerance!(settings::QPALMsettings,tolerance::Float)::Nothing
	settings.eps_dual_inf = tolerance
	return
end

# wrapper to set the iteration limit of each solver in a consistent way
function set_iterationsLimit!(settings::QPALMsettings,limit::Int)::Nothing
	# set the appropriate setting to the value
	settings.max_iter=limit
	return
end

# wrapper to set the time limit of each solver in a consistent way
function set_timeLimit!(settings::QPALMsettings,limit::Float)::Nothing
	# set the appropriate setting to the value
	settings.time_limit=limit
	return
end


## Workspace ##########################################################
# structure used for storing data for QPALM solver
mutable struct QPALMworkspace <: SubsolverWorkspace
    # objective
    problem::Problem
    # memory
    model::QPALM.Model
    settings::QPALMsettings
    # outdated flag
    outdated::Bool
end

## Setup & Update ##########################################################
# this function creates an QPALM.Model representing the given CvxQproblem
function setup(problem::Problem,settings::QPALMsettings;withPrecompilation::Bool=false,nodeType::Type=AbstractNode)::QPALMworkspace

    # check the problem
    @assert problem.objFun isa NullObjective || problem.objFun isa LinearObjective || problem.objFun isa QuadraticObjective
    @assert problem.cnsSet isa NullConstraintSet || problem.cnsSet isa LinearConstraintSet

    # reformat the settings
    settings_dict = Dict{Symbol,Any}()
    for field in fieldnames(QPALMsettings)
        settings_dict[field] = getfield(settings,field)
    end

    # ensure type consistency
    objFun = QuadraticObjective{SpMatrix{Float},Vector{Float}}(problem.objFun)
    cnsSet = LinearConstraintSet{SpMatrix{Float}}(problem.cnsSet)


    # create the QPALMworkspace
    model = QPALM.Model()
    if length(problem.varSet.loBs) > 0
        QPALM.setup!(model;Q=sparse(objFun.Q),q=objFun.L,
                           A=vcat(speye(get_size(problem.varSet)),sparse(cnsSet.A)),
                           bmin=vcat(problem.varSet.loBs,cnsSet.loBs),
                           bmax=vcat(problem.varSet.upBs,cnsSet.upBs),
                           settings_dict...)
    end

    #create the workspace
	workspace = QPALMworkspace(problem,model,settings,false)

	# precompile the main functions to be used according to the workspace created
	if withPrecompilation
		precompile(make_outdated!,(typeof(workspace),))
		precompile(update!,(typeof(workspace),))
		precompile(solve!,(nodeType,typeof(workspace)))
	end


	return workspace
end

# it marks the workspace as outdated
function make_outdated!(workspace::QPALMworkspace)::Nothing
    workspace.outdated = true
    return
end

#
function update!(workspace::QPALMworkspace)::Nothing

    # reformat the settings
    settings_dict = Dict{Symbol,Any}()
    for field in fieldnames(QPALMsettings)
        settings_dict[field] = getfield(workspace.settings,field)
    end

    # ensure type consistency
    objFun = QuadraticObjective{SpMatrix{Float},Vector{Float}}(workspace.problem.objFun)
    cnsSet = LinearConstraintSet{SpMatrix{Float}}(workspace.problem.cnsSet)

    # setup QPALM for the new problem
    QPALM.setup!(workspace.model;Q=sparse(objFun.Q),q=objFun.L,
                 A=vcat(speye(get_size(workspace.problem.varSet)),sparse(cnsSet.A)),
                 bmin=vcat(workspace.problem.varSet.loBs,cnsSet.loBs),
                 bmax=vcat(workspace.problem.varSet.upBs,cnsSet.upBs),
                 settings_dict...)

    # mark the workspace as up to date
    workspace.outdated = false
    return
end

## Solve ##########################################################
function solve!(node::AbstractNode,workspace::QPALMworkspace;objUpperLimit::Float=Inf)::Tuple{Int8,Float}
	# collect info on the problem
	numVars = get_size(workspace.problem.varSet)
	numCnss = get_size(workspace.problem.cnsSet)
	withCuts = nnz(sparse(node.cutSet.A)) > 0

	# check if local cuts are present
    if withCuts # there are some local cuts

        # construct a temporary problem definition (to accomodate the cuts)
        tmpProblem = deepcopy(workspace.problem)
        append!(tmpProblem.cnsSet,LinearConstraintSet(node.cutSet.A[[1],:],[node.cutSet.loBs[1]],[node.cutSet.upBs[1]]))
		update_bounds!(tmpProblem.varSet,loBs=node.varLoBs,upBs=node.varUpBs)
		update_bounds!(tmpProblem.cnsSet,collect(1:length(node.cnsLoBs)),loBs=node.cnsLoBs,upBs=node.cnsUpBs)

        # construct a temporary QPALM model
        tmpWorkspace = setup(tmpProblem,workspace.settings)
		QPALM.warm_start!(tmpWorkspace.model; x_warm_start=node.primal, y_warm_start=vcat(node.bndDual,node.cnsDual,node.cutDual))

        # solve the temporary problem
        sol = QPALM.solve!(tmpWorkspace.model)

    else # there is no local cut
        # update the problem formulation if needed
        if workspace.outdated
            update!(workspace)
        end

        # update bounds in the the qpalm model
        QPALM.update!(workspace.model;bmin=vcat(node.varLoBs,node.cnsLoBs),bmax=vcat(node.varUpBs,node.cnsUpBs))

        # set hotstart info
        if length(node.primal) > 0 && length(node.bndDual) > 0 && length(node.cnsDual) > 0
            QPALM.warm_start!(workspace.model; x_warm_start=node.primal, y_warm_start=vcat(node.bndDual,node.cnsDual))
        end

        # solve problem
        sol = QPALM.solve!(workspace.model)
    end

    # output sol info
    if  sol.info.status_val == 1
        status = 0 # "solved"
        @. node.primal = sol.x
        @. node.bndDual = sol.y[1:numVars]
        @. node.cnsDual = sol.y[numVars+1:numVars+numCnss]
		if withCuts
			@. node.cutDual = sol.y[numVars+numCnss+1:end]
		end
		oldObjLoB = node.objLoB
		update_objBounds!(node,workspace.problem,workspace.settings.eps_prim_inf,10*workspace.settings.eps_dual_inf)
		node.objLoB = max(node.objLoB,oldObjLoB) # avoid problems with lack of accuracy of the dual
	elseif sol.info.status_val == -3
        status = 1 # "infeasible"
		# update the primal info
        @. node.primal = @. min(max(sol.x,node.varLoBs),node.varUpBs)
		# do not update the dual info because
		# the old ones are still good
		# and we cannot trust the new ones
		if withCuts
			@. node.cutDual = sol.y[numVars+numCnss+1:end]
		end
        node.objLoB = node.objUpB = Inf
	elseif sol.info.status_val == 2 && workspace.settings.max_iter < maxQPALMiterations
		status = 1 # solved
		@. node.primal = min(max(sol.x,node.varLoBs),node.varUpBs)
        @. node.bndDual = sol.y[1:numVars]
        @. node.cnsDual = sol.y[numVars+1:numVars+numCnss]
		if withCuts
			@. node.cutDual = sol.y[numVars+numCnss+1:end]
		end
		update_objBounds!(node,workspace.problem,workspace.settings.eps_prim_inf,10*workspace.settings.eps_dual_inf)

    elseif sol.info.status_val in [2,3,4,-6,-2]
        status = 2 # "unreliable"
        @. node.primal = min(max(sol.x,node.varLoBs),node.varUpBs)
        @. node.bndDual = sol.y[1:numVars]
        @. node.cnsDual = sol.y[numVars+1:numVars+numCnss]
		if withCuts
			@. node.cutDual = sol.y[numVars+numCnss+1:end]
		end
		oldObjLoB = node.objLoB
		update_objBounds!(node,workspace.problem,workspace.settings.eps_prim_inf,10*workspace.settings.eps_dual_inf)
		node.objLoB = max(node.objLoB,oldObjLoB) # it is possible that the parent node solution was more successful
		@warn "Inaccuracy in node sol, status: "*string(sol.info.status)*" (code: "*string(status)*")"
    elseif sol.info.status_val in [-7,-10]
        status = 3 # "error"
        @error "Subsover error, status: "*string(sol.info.status)*" (code: "*string(sol.info.status_val)*")"
    else
		status = 3
        @error "Subsolver unknown status: "*string(sol.info.status)*"("*string(sol.info.status_val)*")"
	end
    return (status, sol.info.run_time)
end
