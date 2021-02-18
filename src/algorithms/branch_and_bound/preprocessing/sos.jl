# @Author: Wim Van Roy
# @Date:   Mi Sep 27 14:26:39 CEST 2019
# @Email:  wim.vanroy@kuleuven.be
# @Filename: sos.jl
# @Last modified by:   massimo
# @Last modified time: 2020-04-28T13:59:53+02:00
# @License: LGPL-3.0
# @Copyright: {{copyright}}


# Perform sos1 groups
function preprocess_sos1!(varsToCheck::Vector{Int}, sos1Groups::Vector{Int}, dscIndices::Array{Int, 1},
                          varLoBs::Vector{Float},varUpBs::Vector{Float})::Tuple{Bool,Set{Int}}

    if 0 in varsToCheck
        groupsToCheck = unique(sos1Groups)
    else
        # Get all dscIndices that belong to a sos1Group
        groupsToCheck = unique(sos1Groups[findall(x->x in varsToCheck, dscIndices)])
    end
    groupsToCheck = filter(x->x!=0, groupsToCheck)

    updatedVars = Set{Int}()

    while length(groupsToCheck) > 0
        groupId = pop!(groupsToCheck)
        ids = dscIndices[findall(x->x==groupId, sos1Groups)]

        # For these ids, check which one is
        positiveIds = ids[varLoBs[ids] .> 0]
        negativeIds = ids[varUpBs[ids] .< 0]
        nonzero = unique(vcat(positiveIds, negativeIds))
        if length(nonzero) > 1
            return false, updatedVars
        elseif length(nonzero) == 1
            ids = filter(x->x!=nonzero[1], ids)
            subIds = filter(x->(varLoBs[x] != 0 || varUpBs[x] != 0), ids)
            union!(updatedVars,subIds)
            varLoBs[subIds] .= 0
            varUpBs[subIds] .= 0
        end
    end

    return true, updatedVars
end
