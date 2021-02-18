# @Author: Massimo De Mauri <massimo>
# @Date:   2020-11-23T11:32:27+01:00
# @Email:  massimo.demauri@gmail.com
# @Filename: HBBworkspace.jl
# @Last modified by:   massimo
# @Last modified time: 2021-01-14T16:11:00+01:00
# @License: LGPL-3.0
# @Copyright: {{copyright}}

mutable struct HBBworkspace <: SupersolverWorkspace
    # problem description
    problem::Problem
    # subsolvers Workspaces
    nlpStepWS::HBBnlpStepWorkspace
    mipStepWS::HBBmipStepWorkspace
    # store all the feasible nodes found in decreasing order: the best first
    feasibleNodes::Vector{BBnode}
    # remember which constraints are linear and which are not
    linearCnssIndices::Vector{Int}
    nonLinearCnssIndices::Vector{Int}
    # hybrid branch and bound status
    status::HBBstatus
    # heuristic
    heuristicsWS::HBBheuristicsWorkspace
    # user settings
    settings::HBBsettings
    # already tested discrete assignments
    blackList::BBblackList
    # workspace status
    outdated::Bool
end
