using Clp; const CLP = Clp.ClpCInterface



mutable struct CLPsettings <: AbstractSettings end
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
	presolveType::Int
	solveType::Int
	infeasibleReturn::Int
end

set_algorithm
set_primal_tolerance
set_dual_tolerance
set_objective_offset
set_dual_objective_limit
set_maximum_iterations
set_maximum_seconds
set_log_level
set_presolve_type
set_solve_type
set_infeasible_return




set_col_solution,
set_column_status,
set_row_status





function CLPsettings(;  algorithm::Int=0,
						primalTolerance::Float64=1e-6,
                        dualTolerance::Float64=1e-6,
						objectiveOffset::Float64=0.0,
                        dualObjectiveLimit::Float64=1.7976931348623157e308,
                        maximumIterations::Int=100000,
                        maximumSeconds::Float64=-1.0,
                        logLevel::Int=0,
                        scalingMode::Int=0,
                        perturbationMode::Int=100,
                        presolveType::Int=0, #TODO review the last settings
                        solveType::Int=0,
                        infeasibleReturn::Int=0)::CLPsettings

   return CLPsettings(algorithm,primalTolerance,dualTolerance,objectiveOffset,
   					  dualObjectiveLimit,maximumIterations,maximumSeconds,
					  logLevel,scalingMode,perturbationMode,
                      presolveType,solveType,infeasibleReturn)::CLPsettings
end

mutable struct CLPworkspace <: AbstractWorkspace
   # problem
   problem::Problem
   # memory
   model::CLP.ClpModel
   settings::CLPsettings
   # outdated flag
   outdated::Bool
end


## Setup & Update ##########################################################
# this function creates an Clp model representing the given linear problem
function setup(problem::Problem,settings::CLPsettings;bb_primalTolerance::Float64=Inf,bb_timeLimit=Inf)::OSQPworkspace

   # check the problem
   @assert problem.objFun isa NullObjective || problem.objFun isa LinearObjective
   @assert problem.cnsSet isa NullConstraintSet || problem.cnsSet isa LinearConstraintSet

   # overwrite the Clp settings depending on the branch and bound settings
    settings.eps_prim_inf = min(settings.primalTolerance,bb_primalTolerance*1e-1)
    if bb_timeLimit < Inf
        if settings.timeLimit <= 0.0
            settings.maximumSeconds = bb_timeLimit
        else
            settings.maximumSeconds = min(settings.maximumSeconds,bb_timeLimit)
        end
    end

    # ensure type consistency
    objFun = LinearObjective{Array{Float64,1}}(problem.objFun)
    cnsSet = LinearConstraintSet{SparseMatrixCSC{Float64,Int64}}(problem.cnsSet)

    # create the Clp workspace
    CLP.load_probem(model,cnsSet.A,problem.varSet.loBs,problem.varSet.upBs,objFun.L,cnsSet.loBs,cnsSet.upBs) # insert the settings somehow

    return CLPworkspace(problem,model,settings,false)
end


# it marks the workspace as outdated
function make_outdated!(workspace::CLPworkspace)::Nothing
    workspace.outdated = true
    return
end


#
function update!(workspace::CLPworkspace)::Nothing

    # reformat the settings
    #TODO


    # ensure type consistency
    objFun = LinearObjective{Array{Float64,1}}(workspace.problem.objFun)
    cnsSet = LinearConstraintSet{SparseMatrixCSC{Float64,Int64}}(workspace.problem.cnsSet)

    # re-setup Clp for the new problem
    CLP.load_problem(model,cnsSet.A,problem.varSet.loBs,problem.varSet.upBs,
                     objFun.L,cnsSet.loBs,cnsSet.upBs) # insert the settings somehow

    # mark the workspace as up to date
    workspace.outdated = false
    return
end


## Solve ##########################################################
function solve!(node::BBnode,workspace::CLPworkspace)::Tuple{Int8,Float64}

   # collect info on the problem
	numVars = get_size(workspace.problem.varSet)
	numCnss = get_size(workspace.problem.cnsSet)
	withCuts = nnz(node.cuts) > 0

	#TODO adapt the following lines
	# check if local cuts are present
	if withCuts # there are some local cuts


		# ensure type consistency
	    objFun = LinearObjective{Array{Float64,1}}(problem.objFun)
	    cnsSet = LinearConstraintSet{SparseMatrixCSC{Float64,Int64}}(problem.cnsSet)

	    # construct a CLP model

	    CLP.load_probem(model,cnsSet.A,problem.varSet.loBs,problem.varSet.upBs,objFun.L,cnsSet.loBs,cnsSet.upBs) # insert the settings somehow

		tmpProblem = deepcopy(workspace.problem)
		append!(tmpProblem.cnsSet,node.cuts)
		update_bounds!(tmpProblem.varSet,loBs=node.varLoBs,upBs=node.varUpBs)
		update_bounds!(tmpProblem.cnsSet,collect(1:length(node.cnsLoBs)),loBs=node.cnsLoBs,upBs=node.cnsUpBs)
		tmpModel = setup(tmpProblem,workspace.settings).model

		# set hotstart info
		if length(node.primal) > 0 && length(node.bndDual) > 0 && length(node.cnsDual) > 0 && length(node.cutDual)
			CLP.set_col_solution(tmpModel,node.primal)
			CLP.set_column_status(tmpModel,node.bndDual)
			CLP.set_row_status(tmpModel,vcat(node.cnsDual,node.cutDual))
		end

		# solve the temporary problem
		# solve problem
		if settings.algorithm == "primal"
			CLPstatus = CLP.primal(tmpModel,0)
		elseif settings.algorithm == "dual"
			CLPstatus = CLP.dual(tmpModel,0)
		elseif settings.algorithm == "barrier"
			CLPstatus = CLP.barrier(tmpModel,0)
		else
			@error "CLP : algorithm unavailable"
		end

		# copy the results in the node
		@. node.primal = CLP.primal_column_solution(tmpModel)
		@. node.bndDual = CLP.dual_column_status(tmpModel)
		@. node.cnsDual = CLP.dual_row_status(tmpModel)[1:numCnss]
		@. node.cutDual = CLP.dual_row_status(tmpModel)[numCnss+1:end]

	else # there is no local cut
		# update the problem formulation if needed
		if workspace.outdated
		  	update!(workspace)
		end

		# construct a temporary problem definition (to accomodate the node bounds)
		tmpProblem = deepcopy(workspace.problem)
		update_bounds!(tmpProblem.varSet,loBs=node.varLoBs,upBs=node.varUpBs)
		update_bounds!(tmpProblem.cnsSet,loBs=node.cnsLoBs,upBs=node.cnsUpBs)


		# set hotstart info
		if length(node.primal) > 0 && length(node.bndDual) > 0 && length(node.cnsDual) > 0
			CLP.set_col_solution(tmpModel,node.primal)
			CLP.set_column_status(tmpModel,node.bndDual)
			CLP.set_row_status(tmpModel,node.cnsDual)
		end

		# solve problem
		if settings.algorithm == "primal"
			CLPstatus = CLP.primal(Model,0)
		elseif settings.algorithm == "dual"
			CLPstatus = CLP.dual(Model,0)
		elseif settings.algorithm == "barrier"
		else
			@error "CLP : algorithm unavailable"
		end

		# copy the results in the node
		@. node.primal = CLP.primal_column_solution(Model)
		@. node.bndDual = CLP.dual_column_solution(Model)
		@. node.cnsDual = CLP.dual_row_solution(Model)[1:numCnss]
		if withCuts
			@. node.cutDual = CLP.dual_row_solution(Model)[numCnss+1:end]
		end
		node.objVal = CLP.get_obj_value(Model)
		node.ObjGap = #TODO

	end


	if CLPstatus == 0
		status = 0 # solved




	end


end







# normal flow
model = CLP.ClpModel()
CLP.load_problem(model,sparse(ones(1,2)),[-5.,-5.],[5.,5.],[3.,1.],[0.],[0.])
status = CLP.initial_solve(model,0)
primal = CLP.primal_column_solution(model)
bndDual = CLP.dual_column_solution(model)
cnsDual = CLP.dual_row_solution(model)
objVal = CLP.get_obj_value(model)
dualObjVal = CLP.get_

CLP.get_col_solution(model)
dualSolution = CLP.get_row_price(model)


function set_row_price(model::CLP.ClpModel, input::Vector{Float64})::Nothing
   CLP._jl__check_model(model)
   CLP.@clp_ccall setRowPrice Cvoid (Ptr{Cvoid},Ptr{Float64}) model.p input
   return
end





set_do_dual
set_do_singleton
set_do_doubleton
set_do_tripleton
set_do_tighten
set_do_forcing
set_do_implied_free
set_do_dupcol
set_do_duprow
set_do_singleton_column
set_presolve_actions
set_substitution
