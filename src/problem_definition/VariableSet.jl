# @Author: Massimo De Mauri <massimo>
# @Date:   2019-04-23T15:05:12+02:00
# @Email:  massimo.demauri@gmail.com
# @Filename: VariableSet.jl
# @Last modified by:   massimo
# @Last modified time: 2019-11-22T11:58:30+01:00
# @License: LGPL-3.0
# @Copyright: {{copyright}}


# abstract and null types
abstract type AbstractVariableSet end
struct NullVariableSet <: AbstractVariableSet end



# struct to store the variables data
mutable struct VariableSet <: AbstractVariableSet
    loBs::Array{Float64,1}
    upBs::Array{Float64,1}
    vals::Array{Float64,1}
    dscIndices::Array{Int,1}
    sos1Groups::Array{Int,1} # assume group 0 as no group
    pseudoCosts::Tuple{Array{Float64,2},Array{Int,2}}
end

# named constructor
function VariableSet(;loBs::Array{Float64,1},upBs::Array{Float64,1},vals::Array{Float64,1}=Float64[],
                      dscIndices::Array{Int,1}=Int[],sos1Groups::Array{Int,1}=Int[],pseudoCosts::Tuple{Array{Float64,2},Array{Int,2}}=(NaNs(0,2),Int.(zeros(0,2))))::VariableSet

    # check the correctness of inputs
    if length(sos1Groups) == 0
        sos1Groups = repeat([0],length(dscIndices))
    elseif length(sos1Groups) != length(dscIndices)
        @error "sos1Groups should either be empty or have the same length of dscIndices"
    end

    if size(pseudoCosts[1],1) == 0
        pseudoCosts = (1e-4*ones(length(dscIndices),2),pseudoCosts[2])
    end
    if size(pseudoCosts[2],1) == 0
        pseudoCosts = (pseudoCosts[1],Int.(zeros(length(dscIndices),2)))
    end
    if !(size(pseudoCosts[1],1) == size(pseudoCosts[2],1) == length(dscIndices)) || size(pseudoCosts[1],2) != 2 || size(pseudoCosts[2],2) != 2
        @error "pseudoCosts should either be empty or have size: (length(dscIndices) x 2, length(dscIndices) x 2) "
    end

    if length(vals) == 0
        @assert length(loBs) == length(upBs)
        vals = Array{Float64,1}(undef,length(loBs))
        for i in 1:length(loBs)
            scenario = (loBs[i]>-Inf) + 2*(upBs[i]<Inf)
            if scenario == 3
                vals[i] = .5*(loBs[i] + upBs[i])
            elseif scenario == 2
                vals[i] = upBs[i]
            elseif scenario == 1
                vals[i] = loBs[i]
            else
                vals[i] = 0
            end
        end
    else
        @assert length(loBs) == length(upBs) == length(vals)
    end

    return VariableSet(loBs,upBs,vals,dscIndices,sos1Groups,pseudoCosts)
end

# type conversion
function VariableSet(variableSet::VariableSet)::VariableSet
    return variableSet
end


# copy functions (Fundamental. These are used in Branch and Bound)
import Base.copy
function copy(variableSet::VariableSet)::VariableSet
    return VariableSet(variableSet.loBs,variableSet.upBs,variableSet.vals,
                       variableSet.dscIndices,variableSet.sos1Groups,variableSet.pseudoCosts)
end
import Base.deepcopy
function deepcopy(variableSet::VariableSet)::VariableSet
    return VariableSet(copy(variableSet.loBs),copy(variableSet.upBs),copy(variableSet.vals),
                       copy(variableSet.dscIndices),copy(variableSet.sos1Groups),deepcopy(variableSet.pseudoCosts))
end


# inspect functions (Fundamental. These are used in Branch and Bound)
function get_size(variableSet::VariableSet)::Int
    return length(variableSet.loBs)
end

function get_numDiscrete(variableSet::VariableSet)::Int
    return length(variableSet.dscIndices)
end

function get_bounds(variableSet::VariableSet)::Tuple{Array{Float64,1},Array{Float64,1}}
    return (variableSet.loBs,variableSet.upBs)
end

function get_discreteIndices(variableSet::VariableSet)::Array{Int,1}
    return variableSet.dscIndices
end

function get_sos1Groups(variableSet::VariableSet)::Array{Int,1}
    return variableSet.sos1Groups
end

function get_pseudoCosts(variableSet::VariableSet)::Tuple{Array{Float64,2},Array{Int,2}}
    return variableSet.pseudoCosts
end


# update functions (Not fundamental. These are used only during problem update)
function remove_variables!(variableSet::VariableSet,indices::Array{Int,1})::Nothing
    # collect info
    varToKeep = filter(x->!(x in indices), collect(1:get_size(variableSet)))
    dscToKeep = [i for i in 1:length(variableSet.dscIndices) if !(variableSet.dscIndices[i] in indices)]
    dscMask = Array{Bool,1}(undef,get_size(variableSet));
    @. dscMask = false; @. dscMask[variableSet.dscIndices] = true
    # eliminate the variables
    variableSet.loBs = variableSet.loBs[varToKeep]
    variableSet.upBs = variableSet.upBs[varToKeep]
    variableSet.vals = variableSet.vals[varToKeep]
    variableSet.dscIndices = findall(dscMask[varToKeep])
    variableSet.sos1Groups = variableSet.sos1Groups[dscToKeep]
    variableSet.pseudoCosts = (variableSet.pseudoCosts[1][dscToKeep,:],
                               variableSet.pseudoCosts[2][dscToKeep,:])

    return
end


function insert_variables!(variableSet1::VariableSet,variableSet2::AbstractVariableSet,insertionPoint::Int)::Nothing
    if variableSet2 isa NullVariableSet
        return
    else
        # collect info
        dscInsertionPoint = findfirst(variableSet1.dscIndices .>=insertionPoint)
        if isnothing(dscInsertionPoint)
            dscInsertionPoint = length(variableSet1.dscIndices)+1
        end
        numNewVariables = length(variableSet2.loBs)
        numNewDiscreteVariables = length(variableSet2.dscIndices)
        sos1Offset = maximum(variableSet1.sos1Groups)

        # append variables
        splice!(variableSet1.loBs,insertionPoint:insertionPoint-1,copy(variableSet2.loBs))
        splice!(variableSet1.upBs,insertionPoint:insertionPoint-1,copy(variableSet2.upBs))
        splice!(variableSet1.vals,insertionPoint:insertionPoint-1,copy(variableSet2.vals))

        splice!(variableSet1.dscIndices,dscInsertionPoint:dscInsertionPoint-1,copy(variableSet2.dscIndices))
        @. variableSet1.dscIndices[dscInsertionPoint:dscInsertionPoint+numNewDiscreteVariables-1] += insertionPoint - 1
        @. variableSet1.dscIndices[dscInsertionPoint+numNewDiscreteVariables:end] += numNewVariables
        splice!(variableSet1.sos1Groups,dscInsertionPoint:dscInsertionPoint-1,@. variableSet2.sos1Groups + sos1Offset*(variableSet2.sos1Groups!=0))
        variableSet1.pseudoCosts = (vcat(variableSet1.pseudoCosts[1][1:dscInsertionPoint-1,:],variableSet2.pseudoCosts[1],variableSet1.pseudoCosts[1][dscInsertionPoint:end,:]),
                                    vcat(variableSet1.pseudoCosts[2][1:dscInsertionPoint-1,:],variableSet2.pseudoCosts[2],variableSet1.pseudoCosts[2][dscInsertionPoint:end,:]))
        return
    end
end


function append_variables!(variableSet1::VariableSet,variableSet2::AbstractVariableSet)::Nothing
    if variableSet2 isa NullVariableSet
        return
    else
        # collect info
        indicesOffset = get_size(variableSet1)
        sos1Offset = maximum(variableSet1.sos1Groups)
        # append variables
        variableSet1.loBs = vcat(variableSet1.loBs,variableSet2.loBs)
        variableSet1.upBs = vcat(variableSet1.upBs,variableSet2.upBs)
        variableSet1.vals = vcat(variableSet1.vals,variableSet2.vals)
        variableSet1.dscIndices = vcat(variableSet1.dscIndices,@. variableSet2.dscIndices + indicesOffset)
        variableSet1.sos1Groups = vcat(variableSet1.sos1Groups,@. variableSet2.sos1Groups + sos1Offset*(variableSet2.sos1Groups!=0))
        variableSet1.pseudoCosts = (vcat(variableSet1.pseudoCosts[1],variableSet2.pseudoCosts[1]),
                                    vcat(variableSet1.pseudoCosts[2],variableSet2.pseudoCosts[2]))
        return
    end
end

function update_bounds!(variableSet::VariableSet;loBs::Array{Float64,1}=Float64[],upBs::Array{Float64,1}=Float64)::Nothing
    if length(loBs) > 0
        @assert length(loBs) == length(variableSet.loBs) == length(variableSet.upBs)
        @. variableSet.loBs = loBs
    end
    if length(upBs) > 0
        @assert length(upBs) == length(variableSet.loBs) == length(variableSet.upBs)
        @. variableSet.upBs = upBs
    end
    return
end

function update_bounds!(variableSet::VariableSet,indices::Array{Int,1};loBs::Array{Float64,1}=Float64[],upBs::Array{Float64,1}=Float64)::Nothing
    if length(loBs) > 0
        @assert length(loBs) == length(indices)
        @. variableSet.loBs[indices] = loBs
    end
    if length(upBs) > 0
        @assert length(upBs) == length(indices)
        @. variableSet.upBs[indices] = upBs
    end
    return
end
