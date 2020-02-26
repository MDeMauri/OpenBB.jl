# @Author: Massimo De Mauri <massimo>
# @Date:   2019-02-19T20:29:29+01:00
# @Email:  massimo.demauri@gmail.com
# @Filename: OSQP_interface.jl
# @Last modified by:   massimo
# @Last modified time: 2020-02-19T14:41:57+01:00
# @License: LGPL-3.0
# @Copyright: {{copyright}}

using OSQP

## Settings ##########################################################
# structure used to hold the settings for OSQP
mutable struct OSQPsettings <: AbstractSettings
    rho::Float64                        # ADMM rho step 	0 < rho
    sigma::Float64                      # ADMM sigma step 	0 < sigma
    max_iter::Int                       # Maximum number of iterations 	0 < max_iter (integer)
    eps_abs::Float64                    # Absolute tolerance 	0 <= eps_abs
    eps_rel::Float64                    # Relative tolerance 	0 <= eps_rel
    eps_prim_inf::Float64               # Primal infeasibility tolerance 	0 <= eps_prim_inf
    eps_dual_inf::Float64               # Dual infeasibility tolerance 	0 <= eps_dual_inf
    alpha::Float64                      # ADMM overrelaxation parameter 	0 < alpha < 2
    delta::Float64                      # Polishing regularization parameter 	0 < delta
    polish::Bool                        # Perform polishing 	True/False
    polish_refine_iter::Int             # Refinement iterations in polish 	0 < polish_refine_iter
    verbose::Bool                       # Print output 	True/False
    scaled_termination::Bool            # Scaled termination conditions 	True/False
    check_termination::Int              # Check termination interval 	0 (disabled) or 0 < check_termination
    warm_start::Bool                    # Perform warm starting 	True/False
    scaling::Int                        # Number of scaling iterations 	0 (disabled) or 0 < scaling (integer)
    adaptive_rho::Bool                  # Adaptive rho 	True/False
    adaptive_rho_interval::Int          # Adaptive rho interval 	0 (automatic) or 0 < adaptive_rho_interval
    adaptive_rho_tolerance::Float64     # Tolerance for adapting rho 	1 <= adaptive_rho_tolerance
    adaptive_rho_fraction::Float64      # Adaptive rho interval as fraction of setup time (auto mode) 	0 < adaptive_rho_fraction
    timeLimit::Float64	                # Run time limit in seconds 	0 (disabled) or 0 <= timeLimit

end


function OSQPsettings(; rho::Float64=1e-1,
                        sigma::Float64=1e-6,
                        max_iter::Int=10000,
                        eps_abs::Float64=1e-6,
                        eps_rel::Float64=1e-6,
                        eps_prim_inf::Float64=1e-4,
                        eps_dual_inf::Float64=1e-4,
                        alpha::Float64=1.6,
                        delta::Float64=1e-06,
                        polish::Bool=true,
                        polish_refine_iter::Int=15,
                        verbose::Bool=false,
                        scaled_termination::Bool=false,
                        check_termination::Int=15,
                        warm_start::Bool=true,
                        scaling::Int=15,
                        adaptive_rho::Bool=true,
                        adaptive_rho_interval::Int=0,
                        adaptive_rho_tolerance::Float64=5.,
                        adaptive_rho_fraction::Float64=0.4,
                        timeLimit::Float64=0.)::OSQPsettings



    return OSQPsettings(rho,sigma,max_iter,eps_abs,eps_rel,eps_prim_inf,eps_dual_inf,alpha,delta,
                        polish,polish_refine_iter,verbose,scaled_termination,check_termination,
                        warm_start,scaling,adaptive_rho,adaptive_rho_interval,adaptive_rho_tolerance,
                        adaptive_rho_fraction,timeLimit)
end




## Workspace ##########################################################
# structure used for storing data for OSQP solver
mutable struct OSQPworkspace <: AbstractWorkspace
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
function setup(problem::Problem,settings::OSQPsettings;bb_primalTolerance::Float64=Inf,bb_timeLimit=Inf)::OSQPworkspace

    # check the problem
    @assert problem.objFun isa NullObjective || problem.objFun isa LinearObjective || problem.objFun isa QuadraticObjective
    @assert problem.cnsSet isa NullConstraintSet || problem.cnsSet isa LinearConstraintSet

    # overwrite the osqp settings depending on the branch and bound settings
    settings.eps_prim_inf = min(settings.eps_prim_inf,bb_primalTolerance*1e-1)
    if bb_timeLimit < Inf
        if settings.timeLimit == 0.
            settings.timeLimit = bb_timeLimit
        else
            settings.timeLimit = min(settings.timeLimit,bb_timeLimit)
        end
    end

    # reformat the settings
    settings_dict = Dict{Symbol,Any}()
    for field in fieldnames(OSQPsettings)
        settings_dict[field] = getfield(settings,field)
    end

    # ensure type consistency
    objFun = QuadraticObjective{SparseMatrixCSC{Float64,Int64},Array{Float64,1}}(problem.objFun)
    cnsSet = LinearConstraintSet{SparseMatrixCSC{Float64,Int64}}(problem.cnsSet)

    # create the OSQPworkspaces
    model = OSQP.Model()
    OSQP.setup!(model;P=sparse(objFun.Q),q=objFun.L,
                      A=vcat(speye(get_size(problem.varSet)),sparse(cnsSet.A)),
                      l=vcat(problem.varSet.loBs,cnsSet.loBs),
                      u=vcat(problem.varSet.upBs,cnsSet.upBs),
                      settings_dict...)


    return OSQPworkspace(problem,model,settings,false)
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
    objFun = QuadraticObjective{SparseMatrixCSC{Float64,Int64},Array{Float64,1}}(workspace.problem.objFun)
    cnsSet = LinearConstraintSet{SparseMatrixCSC{Float64,Int64}}(workspace.problem.cnsSet)

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
function solve!(node::BBnode,workspace::OSQPworkspace)::Tuple{Int8,Float64}

	# collect info on the problem
	numVars = get_size(workspace.problem.varSet)
	numCnss = get_size(workspace.problem.cnsSet)
	withCuts = nnz(sparse(node.cuts.A)) > 0

    # check if local cuts are present
    if withCuts # there are some local cuts

        # construct a temporary problem definition (to accomodate the cuts)
        tmpProblem = deepcopy(workspace.problem)
        append!(tmpProblem.cnsSet,node.cuts)
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
        node.objVal = sol.info.obj_val
        node.objGap = max(workspace.settings.eps_abs,
                          workspace.settings.eps_rel*abs(node.objVal))
    elseif sol.info.status_val == -3
        status = 1 # "infeasible"
        @. node.primal = @. min(max(sol.x,node.varLoBs),node.varUpBs)
        @. node.bndDual = sol.y[1:numVars]
        @. node.cnsDual = sol.y[numVars+1:numVars+numCnss]
		if withCuts
			@. node.cutDual = sol.y[numVars+numCnss+1:end]
		end
        node.objVal = Inf
        node.objGap = 0.0
    elseif sol.info.status_val in [2,4,3,-6,-2]
        status = 2 # "unreliable
        @. node.primal = min(max(sol.x,node.varLoBs-workspace.settings.eps_prim_inf),node.varUpBs+workspace.settings.eps_prim_inf)
        @. node.bndDual = sol.y[1:numVars]
        @. node.cnsDual = sol.y[numVars+1:numVars+numCnss]
		if withCuts
			@. node.cutDual = sol.y[numVars+numCnss+1:end]
		end
        objFun = QuadraticObjective{SparseMatrixCSC{Float64,Int64},Array{Float64,1}}(workspace.problem.objFun)
        newObjVal = 1/2 * transpose(node.primal) * objFun.Q *node.primal + transpose(objFun.L) * node.primal
        if newObjVal >= node.objVal - node.objGap
            node.objGap = newObjVal - node.objVal + node.objGap #TODO: recopute the gap if possible
            node.objVal = newObjVal
        else
            node.objGap = Inf #TODO: recopute the gap if possible
            @warn "Inaccuracy in node sol, status: "*string(sol.info.status)*" (code: "*string(status)*")"
        end


        @warn "Inaccuracy in node sol, message: "*string(sol.info.status)*" (code: "*string(sol.info.status_val)*")"
    elseif sol.info.status_val in [-7,-10]
        status = 3 # "error"
        @error "Subsover error, status: "*string(sol.info.status)*" (code: "*string(sol.info.status_val)*")"
    else
        @error "Subsolver unknown status: "*string(sol.info.status)*"("*string(sol.info.status_val)*")"
    end

    return (status, sol.info.run_time)
end
