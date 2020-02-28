# @Author: Massimo De Mauri <massimo>
# @Date:   2019-03-12T15:41:14+01:00
# @Email:  massimo.demauri@gmail.com
# @Filename: simple_rounding_heuristics.jl
# @Last modified by:   massimo
# @Last modified time: 2020-02-28T14:46:29+01:00
# @License: LGPL-3.0
# @Copyright: {{copyright}}


function simple_rounding_heuristics(node::BBnode,workspace::BBworkspace)::Array{BBnode,1}

    # round the primal info and fix the discrete variables
    hNode = deepcopy(node)
    hNode.primal[workspace.problem.varSet.dscIndices]  =
    hNode.varLoBs[workspace.problem.varSet.dscIndices] =
    hNode.varUpBs[workspace.problem.varSet.dscIndices] = round.(hNode.primal[workspace.problem.varSet.dscIndices])

    if  all(@. workspace.problem.varSet.loBs <= hNode.primal <= workspace.problem.varSet.upBs)
        # return the resulting node
        return [hNode]
    else
        return BBnode[]
    end
end
