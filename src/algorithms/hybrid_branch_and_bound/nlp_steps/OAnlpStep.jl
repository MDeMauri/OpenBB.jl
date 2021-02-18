# @Author: Massimo De Mauri <massimo>
# @Date:   2020-11-24T11:24:11+01:00
# @Email:  massimo.demauri@gmail.com
# @Filename: classical_nlpStep.jl
# @Last modified by:   massimo
# @Last modified time: 2021-02-16T13:41:33+01:00
# @License: LGPL-3.0
# @Copyright: {{copyright}}

mutable struct OAnlpStepSettings <: HBBnlpStepSettings
    nlpSettings::AbstractSettings
end

# named constructor
function OAnlpStepSettings(;nlpSettings::AbstractSettings=IPOPTsettings())::OAnlpStepSettings
    return OAnlpStepSettings(nlpSettings)
end


mutable struct OAnlpStepWorkspace <: HBBnlpStepWorkspace
    # reference to the general problem (check for changes)
    problem::Problem
    # info for the nlpSubsolver
    nlpSolverWS::Vector{SubsolverWorkspace}
    # processes to be used
    numProcesses::Int
	# store the settings
	settings
    # tells if the workspace needs an update
    outdated::Bool
end

function setup(problem::Problem,stepSettings::OAnlpStepSettings,numProcesses::Int=1;localOnly::Bool=false)::OAnlpStepWorkspace

    @assert numProcesses > 0 && numProcesses <= nprocs()

    if numProcesses > 1 && !localOnly # multi-process setup
        @sync for id in workers()[1:numProcesses-1]
            @async remotecall_fetch(OpenBB.eval,id,:(nlpStepWS = setup($problem,$stepSettings,$numProcesses,localOnly=true)))
        end
    end
    problem_ = deepcopy(problem)
    nlpSolverWS = [setup(problem_,stepSettings.nlpSettings),setup(build_restorationProblem(problem_),stepSettings.nlpSettings)]
    return OAnlpStepWorkspace(problem,nlpSolverWS,numProcesses,stepSettings,false)
end


function update!(workspace::OAnlpStepWorkspace)::Nothing

    @assert workspace.numProcesses > 0 && workspace.numProcesses <= nprocs()

    if workspace.numProcesses > 1 # multi-process setup
        @sync for id in workers()[1:workspace.numProcesses-1]
            @async remotecall_fetch(OpenBB.eval,id,:(nlpStepWS = setup($workspace.problem,$workspace.settings,1)))
        end
    end

    problem_ = deepcopy(workspace.problem)
    workspace.nlpSolverWS = [setup(problem_,workspace.settings.nlpSettings),
    						 setup(build_restorationProblem(problem_),workspace.settings.nlpSettings)]

    return
end


function make_outdated!(workspace::OAnlpStepWorkspace)::Nothing
    workspace.outdated = true
    return
end





function solve!(nodes::VirtualVector{BBnode},workspace::OAnlpStepWorkspace)::Float

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


function local_solve!(nodes::VirtualVector{BBnode},workspace::OAnlpStepWorkspace)::Float

	results = Vector{Tuple{Int8,Float}}(undef,length(nodes))
	for (k,node) in enumerate(nodes)
		# fix the discrete variables of the problem
		dscIndices = workspace.problem.varSet.dscIndices
		@. node.varLoBs[dscIndices] = node.varUpBs[dscIndices] = round(node.primal[dscIndices])
		results[k] = solve!(node,workspace.nlpSolverWS[1])

		if results[k][1] == 2
			println(node.primal[dscIndices])
		end

		# try to restore solution
		if results[k][1] == 2
			@info :restoration
			numCnss = get_size(workspace.problem.cnsSet)
			cnsVals = evaluate(workspace.problem.cnsSet,node.primal)
			violationUpBs = max.(cnsVals-workspace.problem.cnsSet.upBs,0)
			violationLoBs = max.(workspace.problem.cnsSet.loBs-cnsVals,0)
			resNode = BBnode(vcat(node.varLoBs,zeros(2*numCnss)),vcat(node.varLoBs,Infs(2*numCnss)),
							 node.cnsLoBs,node.cnsUpBs,
							 vcat(node.primal,violationUpBs,violationLoBs),
							 vcat(node.bndDual,zeros(2*numCnss)),node.cnsDual)
			results[k] = solve!(resNode,workspace.nlpSolverWS[2])
			node.primal .= resNode.primal[1:length(node.primal)]
			update_objBounds!(node,workspace.problem,get_primalTolerance(workspace.nlpSolverWS[1].settings),get_dualTolerance(workspace.nlpSolverWS[1].settings))
			@info node.objUpB
		end
	end

	return minimum(x->x.objUpB,nodes)
end

function remote_solve!(nodes::Vector{BBnode},workspace::OAnlpStepWorkspace)::Tuple{Float,Vector{BBnode}}
	return (local_solve!(nodes,workspace),nodes)
end




# helper function
function build_restorationProblem(masterProblem::Problem{<:ObjectiveFunction,ConvexConstraintSet{Th,Tj}})::Problem where {Th<:Union{SpMatrix{Float},Matrix{Float}},Tj<:Union{SpMatrix{Float},Matrix{Float}}}
    # construct the infeasibility reduction problem
    # 		min 	w'*s
    # 		s.t. 	g(x,y) <= s
    # where w is a set of weights, but ipopt does it automatically
    varSet,cnsSet = deepcopy(masterProblem.varSet),deepcopy(masterProblem.cnsSet)
    numVars,numCnss = get_size(varSet),get_size(cnsSet)

    # create a set of slack variables
    append!(varSet,VariableSet(vals=zeros(2*numCnss),loBs=zeros(2*numCnss),upBs=Infs(2*numCnss)))

    # relax the constraints using the slack variables
    oldEvalVal = cnsSet.evalVal
    cnsSet.evalVal = function newEvalVal(x::VirtualVector{Float})::Vector{Float}
    					return oldEvalVal(view(x,1:numVars)) - view(x,numVars+1:numVars+numCnss) + view(x,numVars+numCnss+1:numVars+2*numCnss)
    				 end


    oldEvalJcb = cnsSet.evalJcb
    if Tj isa Matrix
        cnsSet.evalJcb = function newEvalJcbD(x::VirtualVector{Float})::Tj
        					return hcat(oldEvalJcb(view(x,1:numVars)),-diagm(ones(numCnss)),diagm(ones(numCnss)))
        				 end
    else
        cnsSet.evalJcb = function newEvalJcbS(x::VirtualVector{Float})::Tj
        					return hcat(oldEvalJcb(view(x,1:numVars)),-spdiagm(0=>ones(numCnss)),spdiagm(0=>ones(numCnss)))
        				 end
    end
    cnsSet.jcbSparsity=hcat(cnsSet.jcbSparsity,speye(Bool,numCnss),speye(Bool,numCnss))

    oldEvalHes = cnsSet.evalHes
    if Th isa Matrix
        cnsSet.evalHes = function newEvalHesD(x::VirtualVector{Float})::Vector{Th}
                            hes_ = oldEvalHes(view(x,1:numVars))
                            for k in 1:numCnss
                                hes_[k] = vcat(hcat(hes_[k],zeros(numVars,2*numCnss)),zeros(2*numCnss,numVars+2*numCnss))
                            end
                            return hes_
                         end
    else
        cnsSet.evalHes = function newEvalHesS(x::VirtualVector{Float})::Vector{Th}
                            hes_ = oldEvalHes(view(x,1:numVars))
                            for k in 1:numCnss
                                hes_[k] = vcat(hcat(hes_[k],spzeros(numVars,2*numCnss)),spzeros(2*numCnss,numVars+2*numCnss))
                            end
                            return hes_
                         end
    end
    for k in 1:numCnss
        cnsSet.hesSparsity[k] = vcat(hcat(cnsSet.hesSparsity[k],spzeros(Bool,numVars,2*numCnss)),spzeros(Bool,2*numCnss,numVars+2*numCnss))
    end

    if Tj isa Matrix
        objFun = LinearObjective(L=vcat(zeros(numVars),ones(2*numCnss)),c=0.0)
    else
        objFun = LinearObjective(L=vcat(spzeros(numVars),spones(2*numCnss)),c=0.0)
    end

    return Problem(varSet=varSet,cnsSet=cnsSet,objFun=objFun)
end

## inspect functions
function get_nlpSettings(workspace::OAnlpStepWorkspace)::SubsolverSettings
	return workspace.nlpSolverWS[1].settings
end
