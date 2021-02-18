# @Author: Massimo De Mauri <massimo>
# @Date:   2019-03-12T15:41:14+01:00
# @Email:  massimo.demauri@gmail.com
# @Filename: simple_rounding_heuristics.jl
# @Last modified by:   massimo
# @Last modified time: 2020-12-17T19:15:40+01:00
# @License: LGPL-3.0
# @Copyright: {{copyright}}


mutable struct BBroundingHeuristicsSettings <: BBheuristicsSettings end

mutable struct BBroundingHeuristicsWS <: BBheuristicsWorkspace
    problem::Problem
    outdated::Bool
end

function setup(problem::Problem,settings::BBroundingHeuristicsSettings,solverSettings::T)::BBroundingHeuristicsWS where T <: AbstractSettings
    return BBroundingHeuristicsWS(problem,false)
end

#
function make_outdated!(heuristicsWS::BBroundingHeuristicsWS)::Nothing
    # mark the heuristic workspace as outdated
    heuristicsWS.outdated = true
    return
end

function update!(heuristicsWS::BBroundingHeuristicsWS)::Nothing
    # mark the heuristic workspace as up-to-date
    heuristicsWS.outdated = false
    return
end


function get_heuristic_nodes(node::BBnode,heuristicsWS::BBroundingHeuristicsWS,workspace::BBworkspace)::Vector{BBnode}

    # update the heuristics heuristicsWS is necessary
    if heuristicsWS.outdated
        update!(heuristicsWS)
    end

    # collect info on the problem
    dscIndices = heuristicsWS.problem.varSet.dscIndices

    # build the heuristic node
    hNode = deepcopy(node)
    @. hNode.primal[dscIndices] = round(node.primal[dscIndices])
    hNode.heuristic = true

    if all(hNode.cnsLoBs .<= evaluate(workspace.problem.cnsSet,hNode.primal) .<= hNode.cnsUpBs)
        # fix the obtained assignment
        hNode.varLoBs[dscIndices] = hNode.varUpBs[dscIndices] = hNode.primal[dscIndices]
        # return the obtained node
        return [hNode]
    else
        # return an empty list (failure)
        return BBnode[]
    end
end
