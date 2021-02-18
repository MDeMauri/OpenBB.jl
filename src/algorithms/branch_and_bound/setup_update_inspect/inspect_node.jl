# @Author: Massimo De Mauri <massimo>
# @Date:   2021-02-12T14:40:54+01:00
# @Email:  massimo.demauri@protonmail.com
# @Filename: inspect_node.jl
# @Last modified by:   massimo
# @Last modified time: 2021-02-12T14:41:00+01:00
# @License: LGPL-3.0
# @Copyright: {{copyright}}

# checks a portion of a solution for feasibility with respect to the given node
function is_feasible(node::BBnode,varsIndices::Array{Int64,1},varsVal::Array{Float64,1},tolerance::Float64)::Bool where T1<:Problem where T2<:AbstractWorkspace where T3<:AbstractSharedMemory
    for (k,index) in enumerate(varsIndices)
        if !(-tolerance + node.varLoBs[index] < varsVal[k] < node.varUpBs[index] + tolerance)
            return false
        end
    end
	return true
end
