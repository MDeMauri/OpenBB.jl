# @Author: Massimo De Mauri <massimo>
# @Date:   2020-01-08T14:39:40+01:00
# @Email:  massimo.demauri@gmail.com
# @Filename: CLP_interface.jl
# @Last modified by:   massimo
# @Last modified time: 2020-01-08T19:03:01+01:00
# @License: LGPL-3.0
# @Copyright: {{copyright}}


using Clp; const CLP = Clp.ClpCInterface



mutable struct CLPsettings <: AbstractSettings
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

function CLPsettings(;  algorithm::Int=0,
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

mutable struct CLPworkspace <: AbstractWorkspace
   # problem
   problem::Problem
   # memory
   settings::CLPsettings
   # outdated flag
   outdated::Bool
end


## Setup & Update ##########################################################
# this function creates an Clp model representing the given linear problem
function setup(problem::Problem,settings::CLPsettings;bb_primalTolerance::Float64=Inf,bb_timeLimit=Inf)::CLPworkspace

   # check the problem
   @assert problem.objFun isa NullObjective || problem.objFun isa LinearObjective
   @assert problem.cnsSet isa NullConstraintSet || problem.cnsSet isa LinearConstraintSet

   # overwrite the Clp settings depending on the branch and bound settings
    settings.primalTolerance = min(settings.primalTolerance,bb_primalTolerance*1e-1)
    if bb_timeLimit < Inf
        if settings.timeLimit <= 0.0
            settings.maximumSeconds = bb_timeLimit
        else
            settings.maximumSeconds = min(settings.maximumSeconds,bb_timeLimit)
        end
    end

    # create the Clp workspace
    return CLPworkspace(problem,settings,false)
end


# it marks the workspace as outdated
function make_outdated!(workspace::CLPworkspace)::Nothing
    workspace.outdated = true
    return
end


#
function update!(workspace::CLPworkspace)::Nothing
    # mark the workspace as up to date
    workspace.outdated = false
    return
end


## Solve ##########################################################
function solve!(node::BBnode,workspace::CLPworkspace)::Tuple{Int8,Float64}

   # collect info on the problem
	numVars = get_size(workspace.problem.varSet)
	numCnss = get_size(workspace.problem.cnsSet)
	withCuts = nnz(sparse(node.cuts.A)) > 0

	# solve the current problem
	tmpModel = CLP.ClpModel()
	CLPstatus = 0

	# ensure type consistency
	objFun = LinearObjective{Array{Float64,1}}(workspace.problem.objFun)
	cnsSet = LinearConstraintSet{SparseMatrixCSC{Float64,Int64}}(workspace.problem.cnsSet)

	# construct a CLP model
	if withCuts # there are some local cuts
		CLP.load_problem(tmpModel,vcat(cnsSet.A,node.cuts.A),
						 node.varLoBs,node.varUpBs,objFun.L,
						 vcat(node.cnsLoBs,node.cuts.loBs),vcat(node.cnsUpBs,node.cuts.upBs))
	else # there is no local cut
		CLP.load_problem(tmpModel,cnsSet.A,
						 node.varLoBs,node.varUpBs,objFun.L,
						 node.cnsLoBs,node.cnsUpBs)
	end

	# set hotstart info
	if length(node.primal) > 0 && length(node.bndDual) > 0 && length(node.cnsDual) > 0
		CLP.set_col_solution(tmpModel,node.primal)
		@. CLP.set_column_status([tmpModel],1:numVars,Int(node.bndDual))
		@. CLP.set_row_status([tmpModel],1:numCnss,Int(node.cnsDual))
	end
	if withCuts && length(node.cutDual) > 0
		@. CLP.set_row_status([tmpModel],numCnss+1:numCnss+get_size(node.cuts),Int(node.cutDual))
	end

	# solve problem
	update_settings!(tmpModel,workspace.settings)
	if workspace.settings.algorithm == 0 # primal
		CLPstatus = CLP.primal(tmpModel,0)
	elseif workspace.settings.algorithm == 1 # dual
		CLPstatus = CLP.dual(tmpModel,0)
	elseif workspace.settings.algorithm == 3 # barrier
		CLPstatus = CLP.initial_barrier_solve(tmpModel,0)
	else
		@error "CLP : algorithm unavailable"
	end


	# copy results in the node
	primal = CLP.primal_column_solution(tmpModel)[1:numVars]
	@. node.primal = primal
	@. node.bndDual = CLP.get_column_status([tmpModel],1:numVars)
	@. node.cnsDual = CLP.get_row_status([tmpModel],1:numCnss)
	if withCuts
		@. node.cutDual =  CLP.get_row_status([tmpModel],numCnss+1:numCnss+get_size(node.cuts))
	end


	if CLPstatus == 0
		status = 0 # solved
		node.objVal = CLP.get_obj_value(tmpModel)
        node.objGap = CLP.dual_tolerance(tmpModel)

	elseif CLPstatus in [1,2]
       status = 1 # infeasible or unbounded
	   node.objVal = Inf
	   node.objGap = 0.0

   	elseif CLPstatus == 3
		status = 2 # unreliable
		node.objVal = CLP.get_obj_value(tmpModel)
		node.objGap = CLP.dual_tolerance(tmpModel) #TODO actually compute the dual bound
		@warn "Inaccuracy in node sol, status code: "*string(CLPstatus)

   	elseif CLPstatus == 4
       status = 3 # "error"
        @error "Subsover error, status code: "*string(CLPstatus)
	end

    return (status,0.0)
end


# updates the CLP settings
function update_settings!(model::CLP.ClpModel,settings::CLPsettings)::Nothing
	CLP.set_primal_tolerance(model,settings.primalTolerance)
	CLP.set_dual_tolerance(model,settings.dualTolerance)
	CLP.set_objective_offset(model,settings.objectiveOffset)
	CLP.set_dual_objective_limit(model,settings.dualObjectiveLimit)
	CLP.set_maximum_iterations(model,settings.maximumIterations)
	CLP.set_maximum_seconds(model,settings.maximumSeconds)
	CLP.set_log_level(model,settings.logLevel)
	CLP.scaling(model,settings.scalingMode)
	CLP.set_perturbation(model,settings.perturbationMode)
	return
end




# stuff I do not know what to do with
# set_do_dual
# set_do_singleton
# set_do_doubleton
# set_do_tripleton
# set_do_tighten
# set_do_forcing
# set_do_implied_free
# set_do_dupcol
# set_do_duprow
# set_do_singleton_column
# set_presolve_actions
# set_substitution

# # normal flow
# model = CLP.ClpModel()
# CLP.load_problem(model,sparse(ones(1,2)),[-5.,-5.],[5.,5.],[3.,1.],[0.],[0.])
# status = CLP.primal(model,0)
# primal = CLP.primal_column_solution(model)
# bndDual = CLP.dual_column_solution(model)
# cnsDual = CLP.dual_row_solution(model)
# objVal = CLP.get_obj_value(model)
# dualObjVal = CLP.get_
#
# CLP.get_col_solution(model)
# dualSolution = CLP.get_row_price(model)
