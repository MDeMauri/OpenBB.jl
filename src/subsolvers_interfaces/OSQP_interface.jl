# @Author: Massimo De Mauri <massimo>
# @Date:   2019-02-19T20:29:29+01:00
# @Email:  massimo.demauri@gmail.com
# @Filename: OSQP_interface.jl
# @Last modified by:   massimo
# @Last modified time: 2021-01-27T17:20:04+01:00
# @License: LGPL-3.0
# @Copyright: {{copyright}}

using OSQP

## Settings ##########################################################
# structure used to hold the settings for OSQP
mutable struct OSQPsettings <: SubsolverSettings
    rho::Float                        	# ADMM rho step 	0 < rho
    sigma::Float                      	# ADMM sigma step 	0 < sigma
    max_iter::Int                       # Maximum number of iterations 	0 < max_iter (integer)
    eps_abs::Float                    	# Absolute tolerance 	0 <= eps_abs
    eps_rel::Float                    	# Relative tolerance 	0 <= eps_rel
    eps_prim_inf::Float               	# Primal infeasibility tolerance 	0 <= eps_prim_inf
    eps_dual_inf::Float               	# Dual infeasibility tolerance 	0 <= eps_dual_inf
    alpha::Float                      	# ADMM overrelaxation parameter 	0 < alpha < 2
    delta::Float                      	# Polishing regularization parameter 	0 < delta
    polish::Bool                        # Perform polishing 	True/False
    polish_refine_iter::Int             # Refinement iterations in polish 	0 < polish_refine_iter
    verbose::Bool                       # Print output 	True/False
    scaled_termination::Bool            # Scaled termination conditions 	True/False
    check_termination::Int              # Check termination interval 	0 (disabled) or 0 < check_termination
    warm_start::Bool                    # Perform warm starting 	True/False
    scaling::Int                        # Number of scaling iterations 	0 (disabled) or 0 < scaling (integer)
    adaptive_rho::Bool                  # Adaptive rho 	True/False
    adaptive_rho_interval::Int          # Adaptive rho interval 	0 (automatic) or 0 < adaptive_rho_interval
    adaptive_rho_tolerance::Float     	# Tolerance for adapting rho 	1 <= adaptive_rho_tolerance
    adaptive_rho_fraction::Float      	# Adaptive rho interval as fraction of setup time (auto mode) 	0 < adaptive_rho_fraction
    timeLimit::Float	                # Run time limit in seconds 	0 (disabled) or 0 <= timeLimit

end


function OSQPsettings(; rho::Float=1e-1,
                        sigma::Float=1e-6,
                        max_iter::Int=10000,
                        eps_abs::Float=1e-6,
                        eps_rel::Float=1e-6,
                        eps_prim_inf::Float=1e-4,
                        eps_dual_inf::Float=1e-6,
                        alpha::Float=1.6,
                        delta::Float=1e-06,
                        polish::Bool=true,
                        polish_refine_iter::Int=15,
                        verbose::Bool=false,
                        scaled_termination::Bool=false,
                        check_termination::Int=15,
                        warm_start::Bool=true,
                        scaling::Int=15,
                        adaptive_rho::Bool=true,
                        adaptive_rho_interval::Int=0,
                        adaptive_rho_tolerance::Float=5.,
                        adaptive_rho_fraction::Float=0.4,
                        timeLimit::Float=0.)::OSQPsettings



    return OSQPsettings(rho,sigma,max_iter,eps_abs,eps_rel,eps_prim_inf,eps_dual_inf,alpha,delta,
                        polish,polish_refine_iter,verbose,scaled_termination,check_termination,
                        warm_start,scaling,adaptive_rho,adaptive_rho_interval,adaptive_rho_tolerance,
                        adaptive_rho_fraction,timeLimit)
end



## Utility Functions ##########################################################

# wrapper to get verbosity of each solver in the same way
function get_verbosity(settings::OSQPsettings)::Bool
	return settings.verbose
end

# wrapper to get primal tolerance of each solver in the same way
function get_primalTolerance(settings::OSQPsettings)::Float
	return settings.eps_prim_inf
end

# wrapper to get dual tolerance of each solver in the same way
function get_dualTolerance(settings::OSQPsettings)::Float
	return settings.eps_dual_inf
end

# wrapper to get iterations limit of each solver in the same way
function get_iterationsLimit(settings::OSQPsettings)::Int
	return settings.max_iter
end

# wrapper to get time limit of each solver in the same way
function get_timeLimit(settings::OSQPsettings)::Float
	return settings.timeLimit
end

# wrapper to set verbosity of each solver in a consistent way
function set_verbosity!(settings::OSQPsettings,verbose::Bool)::Nothing
	settings.verbose = verbose
	return
end

# wrapper to set primal tolerance of each solve in a consistent way
function set_primalTolerance!(settings::OSQPsettings,tolerance::Float)::Nothing
	settings.eps_prim_inf= tolerance
	return
end

# wrapper to set dual tolerance of each solve in a consistent way
function set_dualTolerance!(settings::OSQPsettings,tolerance::Float)::Nothing
	settings.eps_dual_inf = tolerance
	return
end

# wrapper to set the iteration limit of each solver in a consistent way
function set_iterationsLimit!(settings::OSQPsettings,limit::Int)::Nothing
	# set the appropriate setting to the value
	settings.max_iter=limit
	return
end

# wrapper to set the time limit of each solver in a consistent way
function set_timeLimit!(settings::OSQPsettings,limit::Float)::Nothing
	# set the appropriate setting to the value
	settings.timeLimit=limit
	return
end


## Workspace ##########################################################
# structure used for storing data for OSQP solver
mutable struct OSQPworkspace <: SubsolverWorkspace
    # problem
    problem::Problem
    # memory
    model::OSQP.Model
    settings::OSQPsettings
    # outdated flag
    outdated::Bool
end


## Setup & Update ##########################################################
# this function creates an OSQP.Model representing the given CvxQproblem
function setup(problem::Problem,settings::OSQPsettings;withPrecompilation::Bool=false,nodeType::Type=AbstractNode)::OSQPworkspace

    # check the problem
    @assert problem.objFun isa NullObjective || problem.objFun isa LinearObjective || problem.objFun isa QuadraticObjective
    @assert problem.cnsSet isa NullConstraintSet || problem.cnsSet isa LinearConstraintSet

    # reformat the settings
    settings_dict = Dict{Symbol,Any}()
    for field in fieldnames(OSQPsettings)
        settings_dict[field] = getfield(settings,field)
    end

    # ensure type consistency
    objFun = QuadraticObjective{SpMatrix{Float},Vector{Float}}(problem.objFun)
    cnsSet = LinearConstraintSet{SpMatrix{Float}}(problem.cnsSet)

    # create the OSQPworkspaces
    model = OSQP.Model()
    OSQP.setup!(model;P=sparse(objFun.Q),q=objFun.L,
                      A=vcat(speye(get_size(problem.varSet)),sparse(cnsSet.A)),
                      l=vcat(problem.varSet.loBs,cnsSet.loBs),
                      u=vcat(problem.varSet.upBs,cnsSet.upBs),
                      settings_dict...)

	#create the workspace
    workspace = OSQPworkspace(problem,model,settings,false)

	# precompile the main functions to be used according to the workspace created
	if withPrecompilation
		precompile(make_outdated!,(typeof(workspace),))
		precompile(update!,(typeof(workspace),))
		precompile(solve!,(nodeType,typeof(workspace)))
	end


	return workspace
end

# it marks the workspace as outdated
function make_outdated!(workspace::OSQPworkspace)::Nothing
    workspace.outdated = true
    return
end

#
function update!(workspace::OSQPworkspace)::Nothing

    # reformat the settingss
    settings_dict = Dict{Symbol,Any}()
    for field in fieldnames(OSQPsettings)
        settings_dict[field] = getfield(workspace.settings,field)
    end

    # ensure type consistency
    objFun = QuadraticObjective{SpMatrix{Float},Vector{Float}}(workspace.problem.objFun)
    cnsSet = LinearConstraintSet{SpMatrix{Float}}(workspace.problem.cnsSet)

    # re-setup OSQP for the new problem
    OSQP.setup!(workspace.model;P=sparse(objFun.Q),q=objFun.L,
                A=vcat(speye(get_size(workspace.problem.varSet)),sparse(cnsSet.A)),
                l=vcat(workspace.problem.varSet.loBs,cnsSet.loBs),
                u=vcat(workspace.problem.varSet.upBs,cnsSet.upBs),
                settings_dict...)

    # mark the workspace as up to date
    workspace.outdated = false
    return
end

## Solve ##########################################################
function solve!(node::AbstractNode,workspace::OSQPworkspace;objUpperLimit::Float=Inf)::Tuple{Int8,Float}

	# collect info on the problem
	numVars = get_size(workspace.problem.varSet)
	numCnss = get_size(workspace.problem.cnsSet)
	withCuts = nnz(sparse(node.cutSet.A)) > 0

    # check if local cuts are present
    if withCuts # there are some local cuts

        # construct a temporary problem definition (to accomodate the cuts)
        tmpProblem = deepcopy(workspace.problem)
        append!(tmpProblem.cnsSet,node.cutSet)
		update_bounds!(tmpProblem.varSet,loBs=node.varLoBs,upBs=node.varUpBs)
		update_bounds!(tmpProblem.cnsSet,collect(1:length(node.cnsLoBs)),loBs=node.cnsLoBs,upBs=node.cnsUpBs)

        # construct a temporary OSQP model
        tmpWorkspace = setup(tmpProblem,workspace.settings)
		OSQP.warm_start!(tmpWorkspace.model; x=node.primal, y=vcat(node.bndDual,node.cnsDual,node.cutDual))

        # solve the temporary problem
        sol = OSQP.solve!(tmpWorkspace.model)

    else # there is no local cut
        # update the problem formulation if needed
        if workspace.outdated
            update!(workspace)
        end

        # update the bounds in the osqp model
        OSQP.update!(workspace.model;l=vcat(node.varLoBs,node.cnsLoBs),u=vcat(node.varUpBs,node.cnsUpBs))

        # set hotstart info
        if length(node.primal) > 0 && length(node.bndDual) > 0 && length(node.cnsDual) > 0
            OSQP.warm_start!(workspace.model; x=node.primal, y=vcat(node.bndDual,node.cnsDual))
        end

        # solve problem
        sol = OSQP.solve!(workspace.model)
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
		update_objBounds!(node,workspace.problem,workspace.settings.eps_prim_inf,workspace.settings.eps_dual_inf)
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
    elseif sol.info.status_val in [2,4,3,-6,-2]
        status = 2 # "unreliable
        @. node.primal = min(max(sol.x,node.varLoBs-workspace.settings.eps_prim_inf),node.varUpBs+workspace.settings.eps_prim_inf)
        @. node.bndDual = sol.y[1:numVars]
        @. node.cnsDual = sol.y[numVars+1:numVars+numCnss]
		if withCuts
			@. node.cutDual = sol.y[numVars+numCnss+1:end]
		end
		oldObjLoB = node.objLoB
		update_objBounds!(node,workspace.problem,workspace.settings.eps_prim_inf,workspace.settings.eps_dual_inf)
		node.objLoB = max(node.objLoB,oldObjLoB) # it is possible that the parent node solution was more successful

        @warn "Inaccuracy in node sol, status: "*string(sol.info.status)*" (code: "*string(sol.info.status_val)*")"
    elseif sol.info.status_val in [-7,-10]
        status = 3 # "error"
        @error "Subsover error, status: "*string(sol.info.status)*" (code: "*string(sol.info.status_val)*")"
    else
        @error "Subsolver unknown status: "*string(sol.info.status)*"("*string(sol.info.status_val)*")"
    end

    return (status, sol.info.run_time)
end
