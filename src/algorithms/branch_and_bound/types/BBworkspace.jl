# @Author: Massimo De Mauri <massimo>
# @Date:   2019-05-23T11:26:28+02:00
# @Email:  massimo.demauri@gmail.com
# @Filename: BBworkspace.jl
# @Last modified by:   massimo
# @Last modified time: 2021-02-09T21:33:19+01:00
# @License: LGPL-3.0
# @Copyright: {{copyright}}


# this structure collects all the informations needed to run B&B
mutable struct BBworkspace <: SupersolverWorkspace
    # problem description
    problem::Problem
    # subsolver WS
    subsolverWS::SubsolverWorkspace
    # multiprocessing communication
    sharedMemory::AbstractSharedMemory
    # branch and bound status
    tree::BBtree
    status::BBstatus
    # nodes updates register
    updatesRegister::BBupdatesRegister
    # heuristic
    heuristicsWS::BBheuristicsWorkspace
    # user settings
    settings::BBsettings
    # solutions to be ignored
    blackList::BBblackList
    # workspace status
    outdated::Bool
    infeasiblesToRecover::Bool
    invalidLowerBounds::Bool
end



# construct the root node for a problem of the given dimensions
function BBroot(workspace::BBworkspace)::BBnode
    numVars = get_numVariables(workspace)
    numCnss = get_numConstraints(workspace)
    varBounds = get_variableBounds(workspace)
    cnsBounds = get_constraintBounds(workspace)

    return BBnode(copy(varBounds[1]),copy(varBounds[2]),
                  copy(cnsBounds[1]),copy(cnsBounds[2]),
                  zeros(numVars),zeros(numVars),zeros(numCnss),
                  workspace.settings.maxNumberOfLocalCuts+ is_mixedBinary(workspace.problem.varSet),
                  workspace.updatesRegister.masterVersion)
end

# construct the root node for a problem of the given dimensions (with primal initialization)
function BBroot(workspace::BBworkspace,primal::Vector{Float})::BBnode
    numVars = get_numVariables(workspace)
    numCnss = get_numConstraints(workspace)
    varBounds = get_variableBounds(workspace)
    cnsBounds = get_constraintBounds(workspace)

    return BBnode(copy(varBounds[1]),copy(varBounds[2]),
                  copy(cnsBounds[1]),copy(cnsBounds[2]),
                  copy(primal),zeros(numVars),zeros(numCnss),
                  workspace.settings.maxNumberOfLocalCuts + is_mixedBinary(workspace.problem.varSet),
                  workspace.updatesRegister.masterVersion)
end
