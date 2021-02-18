# @Author: Massimo De Mauri <massimo>
# @Date:   2020-11-24T11:24:11+01:00
# @Email:  massimo.demauri@gmail.com
# @Filename: classical_nlpStep.jl
# @Last modified by:   massimo
# @Last modified time: 2021-01-26T19:26:54+01:00
# @License: LGPL-3.0
# @Copyright: {{copyright}}

mutable struct POAnlpStepSettings <: HBBnlpStepSettings
    nlpSettings::AbstractSettings
end

# named constructor
function POAnlpStepSettings(;nlpSettings::AbstractSettings=IPOPTsettings())::POAnlpStepSettings
    return POAnlpStepSettings(nlpSettings)
end


mutable struct POAnlpStepWorkspace <: HBBnlpStepWorkspace
    # reference to the general problem (check for changes)
    problem::Problem
    # info for the nlpSubsolver
    nlpSolverWS::Vector{SubsolverWorkspace}
	# pointer to problem parameters space
	paramSpaceRef::RefValue{Vector{Float}}
    # processes to be used
    numProcesses::Int
	# store the settings
	settings::POAnlpStepSettings
    # tells if the workspace needs an update
    outdated::Bool
end

function setup(problem::Problem,settings::POAnlpStepSettings,numProcesses::Int=1)::POAnlpStepWorkspace

    @assert numProcesses > 0 && numProcesses <= nprocs()

	out = Ref{POAnlpStepWorkspace}()
	@sync begin
		if numProcesses > 1 && !localOnly # multi-process setup
			for id in workers()[1:numProcesses-1]
				@async remotecall_fetch(OpenBB.eval,id,:(nlpStepWS = local_setup($problem,$settings,$numProcesses)))
			end
		end
		out[] = local_setup(problem,settings,numProcesses)
	end
    return out[]
end

function local_setup(problem::Problem,settings::POAnlpStepSettings,numProcesses::Int=1)::POAnlpStepWorkspace

	# collect info
	numVars = get_size(problem.varSet)
	dscIndices = problem.varSet.dscIndices

	# build proximal objective
	paramSpaceRef = Ref(zeros(length(dscIndices)))
	function prxEvalVal(x::Vector{Float})::Float
		return sum(@. 500.0*(x[dscIndices]-paramSpaceRef[])^2)
	end
	grdSparsity = sparsevec(dscIndices,true,numVars)
	function prxEvalGrd(x::Vector{Float})::SpVector{Float}
		return sparsevec(dscIndices,1000.0*(x[dscIndices].-paramSpaceRef[]),numVars)
	end
	hesSparsity = sparse(dscIndices,dscIndices,true,numVars,numVars)
	function prxEvalHes(x::Vector{Float})::SpMatrix{Float}
		return sparse(dscIndices,dscIndices,1000.0,numVars,numVars)
	end
	prxObjective = ConvexObjective(evalVal=prxEvalVal,evalGrd=prxEvalGrd,evalHes=prxEvalHes,
								   grdSparsity=grdSparsity,hesSparsity=hesSparsity,
								   typeGrd=SpVector{Float},typeHes=SpMatrix{Float})
	prxProblem = Problem(varSet=problem.varSet,cnsSet=problem.cnsSet,objFun=prxObjective)
	prxSettings = deepcopy(settings.nlpSettings)
	set_dualTolerance!(prxSettings,get_primalTolerance(prxSettings)^2)

    nlpSolverWS = [setup(prxProblem,settings.nlpSettings),setup(problem,settings.nlpSettings)]

    return POAnlpStepWorkspace(problem,nlpSolverWS,paramSpaceRef,numProcesses,settings,false)
end

function update!(workspace::POAnlpStepWorkspace)::Nothing

    @assert workspace.numProcesses > 0 && workspace.numProcesses <= nprocs()

	@sync begin
	    if workspace.numProcesses > 1 # multi-process setup
	        for id in workers()[1:workspace.numProcesses-1]
	            @async remotecall_fetch(OpenBB.eval,id,:(nlpStepWS = local_setup($workspace.problem,$workspace.settings,$numProcesses)))
	        end
	    end
	    workspace_ = local_setup(workspace.problem,workspace.settings,workspace.numProcesses)
		workspace.nlpSolverWS = workspace_.nlpSolverWS
		workspace.paramSpaceRef = workspace_.paramSpaceRef
	end
    return
end


function make_outdated!(workspace::POAnlpStepWorkspace)::Nothing
    workspace.outdated = true
    return
end


function solve!(nodes::VirtualVector{BBnode},workspace::POAnlpStepWorkspace)::Float

    # update the workspace if needed
    if workspace.outdated
        update!(workspace)
    end

    # solve each node and collect results
    bestObjVal = Inf
    if workspace.numProcesses>1 && length(nodes)>1

		# solve some nodes on the remote processes
		objVals = Infs(workspace.numProcesses)
		@sync begin

			# solve some nodes in the current process
			if length(nodes) >= workspace.numProcesses
				@async objVals[1] = local_solve!(view(nodes,1:workspace.numProcesses:length(nodes)),workspace)
			end

			# send other nodes to the remote workers
			index = 1
	        for id in workers()[1:workspace.numProcesses-1]
	            if length(nodes) > index
					index += 1
	                @async begin
						 (objVals[index],nodes[index:workspace.numProcesses:end]) =
						 	remotecall_fetch(OpenBB.eval,id,:(remote_solve!($(nodes[index:workspace.numProcesses:end]),nlpStepWS)))
					end
	            end
	        end
		end

		# collect results
	    bestObjVal = minimum(objVals)

    else # solve all nodes locally
		bestObjVal = local_solve!(nodes,workspace)
    end

    return bestObjVal
end


function local_solve!(nodes::VirtualVector{BBnode},workspace::POAnlpStepWorkspace)::Float

	problem = workspace.problem
	settings = workspace.settings
	# collect info
	numVars = get_size(problem.varSet)
	dscIndices = problem.varSet.dscIndices

	results = Vector{Tuple{Int8,Float}}(undef,length(nodes))
	for (k,node) in enumerate(nodes)

		# collect info
		dscIndices = workspace.problem.varSet.dscIndices
		dscValues = round.(node.primal[dscIndices])

		# set the parameters of the proximal problem
		workspace.paramSpaceRef[][1:end] = dscValues

		# solve the proximal problem
		@info (:before,evaluate(workspace.nlpSolverWS[1].problem.objFun,node.primal))
		results[k] = solve!(node,workspace.nlpSolverWS[1])
		@info (:after,evaluate(workspace.nlpSolverWS[1].problem.objFun,node.primal),node.objUpB)

		# check if we got a discrete solution
		feasibleAssignment = true
		for (k,ind) in enumerate(dscIndices)
			if abs(node.primal[ind]-dscValues[k]) > get_primalTolerance(workspace.settings.nlpSettings)
				feasibleAssignment = false
				break
			end
		end

		# if the assignment is feasible, fix the discrete variables to optimize the continuous ones
		if feasibleAssignment
			@info :feasible
			@. node.primal[dscIndices] = node.varLoBs[dscIndices] = node.varUpBs[dscIndices] = dscValues
			@info results[k] = solve!(node,workspace.nlpSolverWS[2])
		else
			results[k] = (1,results[k][2])
			node.objLoB = node.objUpB = Inf
		end
	end

	return minimum(x->x.objUpB,nodes)
end

function remote_solve!(nodes::Vector{BBnode},workspace::POAnlpStepWorkspace)::Tuple{Float,Vector{BBnode}}
	return (local_solve!(nodes,workspace),nodes)
end


## inspect functions
function get_nlpSettings(workspace::POAnlpStepWorkspace)::SubsolverSettings
	return workspace.nlpSolverWS[1].settings
end
