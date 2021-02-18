# @Author: Massimo De Mauri <massimo>
# @Date:   2019-03-12T15:39:53+01:00
# @Email:  massimo.demauri@gmail.com
# @Filename: heuristics.jl
# @Last modified by:   massimo
# @Last modified time: 2021-01-21T10:59:18+01:00
# @License: LGPL-3.0
# @Copyright: {{copyright}}


####################### Functions For No Heuristics Case #################################
# builds a null heuristic that cannot be ever called
function setup(problem::Problem,settings::NullHBBheuristicsSettings,subsolverSettings::AbstractSettings)::NullHBBheuristicsWorkspace
    return NullHBBheuristicsWorkspace(false)
end

# fake outdated marking operation
function make_outdated!(heuristicsWS::NullHBBheuristicsWorkspace)::Nothing
    return
end

# fake update operation
function update!(heuristicsWS::NullHBBheuristicsWorkspace)::Nothing
    return
end

####################### Triggering Conditions #################################
#TODO : give form to solution
# check satisfaction of triggering condition
function check_heuristicsTriggerCondition(solution,workspace::HBBworkspace)::Bool
    return workspace.settings.heuristicsTriggerCondition[1](solution,workspace,workspace.settings.heuristicsTriggerCondition[2]...)
end

# never trigger heuristics
function never_trigger(solution,workspace::HBBworkspace)::Bool
    return false
end
