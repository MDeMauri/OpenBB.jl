# @Author: Massimo De Mauri <massimo>
# @Date:   2019-12-10T22:59:22+01:00
# @Email:  massimo.demauri@gmail.com
# @Filename: feasibility_pump.jl
# @Last modified by:   massimo
# @Last modified time: 2021-02-16T13:40:56+01:00
# @License: LGPL-3.0
# @Copyright: {{copyright}}

# setting type
mutable struct BBfeasibilityPumpSettings <: BBheuristicsSettings
    maxNumberOfIterations::Int
end

# named constructor with defaults
function BBfeasibilityPumpSettings(;maxNumberOfIterations::Int=1)::BBfeasibilityPumpSettings
    return BBfeasibilityPumpSettings(maxNumberOfIterations)
end


# heuristicsWS
mutable struct BBfeasibilityPumpWS{T}<:BBheuristicsWorkspace where T<:SubsolverWorkspace
    problem::Problem
    solverWS::T
    settings::BBfeasibilityPumpSettings
    outdated::Bool
end

#
function setup(problem::Problem,settings::BBfeasibilityPumpSettings,solverSettings::T)::BBfeasibilityPumpWS where T <: AbstractSettings

    # collect info on the problem
    numVars = get_size(problem.varSet)
    numCnss = get_size(problem.cnsSet)
    dscIndices = problem.varSet.dscIndices
    numDscVars = length(dscIndices)

    # make space for the slack variables needed for the L1 norm
    varSet = deepcopy(problem.varSet)
    append!(varSet,VariableSet(loBs=zeros(numDscVars),upBs=Infs(numDscVars),vals=zeros(numDscVars)))

    # extend the constraint set to define the L1 norm
    Atmp = spzeros(numDscVars,numVars); Atmp[:,dscIndices] = speye(numDscVars)
    L1Constraints = LinearConstraintSet(A = vcat(hcat(Atmp,-speye(numDscVars)),hcat(-Atmp,-speye(numDscVars))),
                                        loBs = -Infs(2*numDscVars),upBs = Infs(2*numDscVars))
    cnsSet = deepcopy(problem.cnsSet)
    append_variables!(cnsSet,numDscVars)
    append!(cnsSet,L1Constraints)

    # define the L1 norm objective
    objFun = LinearObjective(L=vcat(spzeros(numVars),ones(numDscVars)))

    # finally construct the L1-distance problem and its solver
    solverWS = setup(Problem(varSet=varSet,cnsSet=cnsSet,objFun=objFun),solverSettings)

    return BBfeasibilityPumpWS(problem,solverWS,settings,false)
end

#
function make_outdated!(heuristicsWS::BBfeasibilityPumpWS)::Nothing
    # mark the heuristic workspace as outdated
    heuristicsWS.outdated = true
    return
end


function update!(heuristicsWS::BBfeasibilityPumpWS)::Nothing

    # collect info on the problem
    numVars = get_size(heuristicsWS.problem.varSet)
    numCnss = get_size(heuristicsWS.problem.cnsSet)
    dscIndices = heuristicsWS.problem.varSet.dscIndices
    numDscVars = length(dscIndices)

    # make space for the slack variables needed for the L1 norm
    varSet = deepcopy(heuristicsWS.problem.varSet)
    append!(varSet,VariableSet(loBs=zeros(numDscVars),upBs=Infs(numDscVars),vals=zeros(numDscVars)))

    # extend the constraint set to define the L1 norm
    Atmp = spzeros(numDscVars,numVars); Atmp[:,dscIndices] = speye(numDscVars)
    L1Constraints = LinearConstraintSet(A = vcat(hcat(Atmp,-speye(numDscVars)),hcat(-Atmp,-speye(numDscVars))),
                                        loBs = -Infs(2*numDscVars),upBs = Infs(2*numDscVars))
    cnsSet = deepcopy(heuristicsWS.problem.cnsSet)
    append_variables!(cnsSet,numDscVars)
    append!(cnsSet,L1Constraints)

    # define the L1 norm objective
    objFun = LinearObjective(L=vcat(spzeros(numVars),ones(numDscVars)))

    # finally regenerate the subsolver workspace
    heuristicsWS.solverWS = setup(Problem(varSet=varSet,cnsSet=cnsSet,objFun=objFun),heuristicsWS.solverWS.settings)

    # mark the heuristic workspace as up to date
    heuristicsWS.outdated = false
    return
end


function get_heuristic_nodes(node::BBnode,heuristicsWS::BBfeasibilityPumpWS,workspace::BBworkspace)::Vector{BBnode}

    # update the heuristicsWS is necessary
    if heuristicsWS.outdated
        update!(heuristicsWS)
    end

    # collect info on the problem
    numVars = get_size(heuristicsWS.problem.varSet)
    numCnss = get_size(heuristicsWS.problem.cnsSet)
    dscIndices = heuristicsWS.problem.varSet.dscIndices
    numDscVars = length(dscIndices)
    primalTolerance = workspace.settings.primalTolerance

    # build the heuristic node
    hNode = deepcopy(node)
    hNode.primal[dscIndices] = round.(node.primal[dscIndices])
    hNode.heuristic = true

    if all(- primalTolerance .+ hNode.cnsLoBs .<= evaluate(workspace.problem.cnsSet,hNode.primal) .<= hNode.cnsUpBs .+ primalTolerance)
        # fix the obtained assignment
        hNode.varLoBs[dscIndices] = hNode.varUpBs[dscIndices] = hNode.primal[dscIndices]
        # return the obtained node
        return [hNode]

    elseif heuristicsWS.settings.maxNumberOfIterations == 0 # no FP iterations allowed

        # return an empty list (failure)
        return BBnode[]

    else

        # extend the heuristic node for the L1-distance problem
        append!(hNode.varLoBs,zeros(numDscVars))
        append!(hNode.varUpBs,Infs(numDscVars))
        append!(hNode.primal,zeros(numDscVars))
        append!(hNode.bndDual,zeros(numDscVars))
        append!(hNode.cnsLoBs,-Infs(2*numDscVars))
        append!(hNode.cnsUpBs,Infs(2*numDscVars))
        append!(hNode.cnsDual,zeros(2*numDscVars))
        append_variables!(hNode.cutSet,numDscVars)

        # allocate memory to check cycling
        lastAssignment = Vector{Float}(undef,numDscVars)

        # iterate
        it = 1
        while true
            # remember the current assignment
            lastAssignment .= hNode.primal[dscIndices]

            # find a new assignment
            hNode.cnsUpBs[numCnss+1:end] = vcat(hNode.primal[dscIndices],-hNode.primal[dscIndices])
            (status,runtime) = solve!(hNode,heuristicsWS.solverWS)
            workspace.status.numExploredNodes += 1

            # if the solver has failed...
            if status != 0
                # return an empty list (failure)
                return BBnode[]
            end

            # propose the rounding of the new solution as heuristic point
            hNode.primal[dscIndices] = round.(hNode.primal[dscIndices])

            # if the new assignment is feasible...
            if all(-primalTolerance .+ hNode.cnsLoBs[1:numCnss] .< evaluate(heuristicsWS.problem.cnsSet,hNode.primal[1:numVars]) .< hNode.cnsUpBs[1:numCnss] .+ primalTolerance)

                # remove the slack variables from the heuristic node
                hNode.varLoBs = hNode.varLoBs[1:numVars]
                hNode.varUpBs = hNode.varUpBs[1:numVars]
                hNode.primal = hNode.primal[1:numVars]
                hNode.bndDual = hNode.bndDual[1:numVars]
                remove_variables!(hNode.cutSet,collect(numVars+1:numVars+numDscVars))

                # remove the L1 constraint from the heuristic node
                hNode.cnsLoBs = hNode.cnsLoBs[1:numCnss]
                hNode.cnsUpBs = hNode.cnsUpBs[1:numCnss]
                hNode.cnsDual = hNode.cnsDual[1:numCnss]

                # fix the obtained feasible assignment
                hNode.varLoBs[dscIndices] = hNode.varUpBs[dscIndices] = hNode.primal[dscIndices]

                # return the obtained node
                return [hNode]

            # if the maximum number of iterations has been reached or the feasibility pump is cycling...
            elseif it == heuristicsWS.settings.maxNumberOfIterations || all(lastAssignment .== hNode.primal[dscIndices])
                # return an empty list (failure)
                return BBnode[]
            else
                # update iterations counter (and iterate again)
                it += 1
            end
        end
    end
end
