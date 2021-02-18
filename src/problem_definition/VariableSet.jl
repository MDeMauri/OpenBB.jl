# @Author: Massimo De Mauri <massimo>
# @Date:   2019-04-23T15:05:12+02:00
# @Email:  massimo.demauri@gmail.com
# @Filename: VariableSet.jl
# @Last modified by:   massimo
# @Last modified time: 2021-02-16T15:36:48+01:00
# @License: LGPL-3.0
# @Copyright: {{copyright}}

# struct to store the variables data
mutable struct VariableSet
    loBs::Vector{Float}
    upBs::Vector{Float}
    vals::Vector{Float}
    dscIndices::Vector{Int}
    bnrIndices::Vector{Int}
    sos1Groups::Vector{Int} # assume group 0 as no group
    pseudoCosts::Tuple{Matrix{Float},Matrix{Int}}
end

# named constructor
function VariableSet(;loBs::Vector{Float},upBs::Vector{Float},vals::Vector{Float}=Float[],
                      dscIndices::Vector{Int}=Int[],sos1Groups::Vector{Int}=Int[],
                      pseudoCosts::Tuple{Matrix{Float},Matrix{Int}}=(NaNs(0,2),Int.(zeros(0,2))))::VariableSet

    @assert length(loBs) == length(upBs)
    numVars = length(loBs)
    numDiscrete = length(dscIndices)

    # check the correctness of inputs
    if length(vals) == 0
        vals = Vector{Float}(undef,numVars)
        for i in 1:numVars
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
        @assert length(vals) == numVars
    end

    if length(sos1Groups) == 0
        sos1Groups = zeros(Int,numDiscrete)
    else
        @assert length(sos1Groups) == numDiscrete
    end



    bnrIndices = Int[]
    binaryDetectorFun = ind -> -1.0 < loBs[ind] && upBs[ind] < 2.0
    bnrIndices = filter(binaryDetectorFun,dscIndices)

    if size(pseudoCosts[1],1) == 0
        pseudoCosts = (1e-4*ones(numDiscrete,2),pseudoCosts[2])
    end
    if size(pseudoCosts[2],1) == 0
        pseudoCosts = (pseudoCosts[1],zeros(Int,numDiscrete,2))
    end
    if !(size(pseudoCosts[1],1) == size(pseudoCosts[2],1) == length(dscIndices)) || size(pseudoCosts[1],2) != 2 || size(pseudoCosts[2],2) != 2
        error("VariableSet: pseudoCosts array should either be empty or have size: (length(dscIndices),2), (length(dscIndices),2)")
    end

    return VariableSet(loBs,upBs,vals,dscIndices,bnrIndices,sos1Groups,pseudoCosts)
end

# empty constructor
function EmptyVariableSet()::VariableSet

    return VariableSet(Vector{Float}(undef,0),Vector{Float}(undef,0),Vector{Float}(undef,0),
                       Vector{Int}(undef,0),Vector{Int}(undef,0),Vector{Int}(undef,0),
                       (zeros(0,2),zeros(Int,0,2)))
end


# type conversion
function VariableSet(varSet::VariableSet)::VariableSet
    return varSet
end


# copy functions (Fundamental. These are used in Branch and Bound)
function Base.copy(varSet::VariableSet)::VariableSet
    return VariableSet(varSet.loBs,varSet.upBs,varSet.vals,
                       varSet.dscIndices,varSet.bnrIndices,
                       varSet.sos1Groups,varSet.pseudoCosts)
end
function Base.deepcopy(varSet::VariableSet)::VariableSet
    return VariableSet(copy(varSet.loBs),copy(varSet.upBs),copy(varSet.vals),
                       copy(varSet.dscIndices),copy(varSet.bnrIndices),
                       copy(varSet.sos1Groups),deepcopy(varSet.pseudoCosts))
end


# inspect functions (Fundamental. These are used in Branch and Bound)
function is_empty(varSet::VariableSet)::Bool
    return length(varSet.loBs) == 0
end

function get_size(varSet::VariableSet)::Int
    return length(varSet.loBs)
end

function get_numDiscrete(varSet::VariableSet)::Int
    return length(varSet.dscIndices)
end

function get_numBinary(varSet::VariableSet)::Int
    return length(varSet.bnrIndices)
end


function get_bounds(varSet::VariableSet)::Tuple{Vector{Float},Vector{Float}}
    return (varSet.loBs,varSet.upBs)
end

function get_discreteIndices(varSet::VariableSet)::Vector{Int}
    return varSet.dscIndices
end

function get_binaryIndices(varSet::VariableSet)::Vector{Int}
    return varSet.bnrIndices
end

function get_sos1Groups(varSet::VariableSet)::Vector{Int}
    return varSet.sos1Groups
end

function get_pseudoCosts(varSet::VariableSet)::Tuple{Matrix{Float},Matrix{Int}}
    return varSet.pseudoCosts
end


# update functions (Not fundamental. These are used only during problem update)
function fix_variables!(varSet::VariableSet,indices::Union{Vector{Int},UnitRange{Int}},values::Vector{Float})::Nothing
    dscIndices_ = findall(ind->ind in varSet.dscIndices,indices)
    @assert maximum(@. abs(values(dscIndices_)-round(values(dscIndices_)))) == 0
    @. varSet.loBs[indices] = varSet.upBs = values
    return
end


function remove_variables!(varSet::VariableSet,indices::Union{Vector{Int},UnitRange{Int}})::Nothing
    # collect info
    varToKeep = filter(x->!(x in indices), collect(1:get_size(varSet)))
    dscToKeep = [i for i in 1:length(varSet.dscIndices) if !(varSet.dscIndices[i] in indices)]
    dscMask = fill(false,get_size(varSet)); @. dscMask[varSet.dscIndices] = true
    bnrMask = fill(false,get_size(varSet)); @. bnrMask[varSet.bnrIndices] = true

    # eliminate the variables
    varSet.loBs = varSet.loBs[varToKeep]
    varSet.upBs = varSet.upBs[varToKeep]
    varSet.vals = varSet.vals[varToKeep]
    varSet.dscIndices = findall(dscMask[varToKeep])
    varSet.bnrIndices = findall(bnrMask[varToKeep])
    varSet.sos1Groups = varSet.sos1Groups[dscToKeep]
    varSet.pseudoCosts = (varSet.pseudoCosts[1][dscToKeep,:],
                          varSet.pseudoCosts[2][dscToKeep,:])

    return
end

function Base.:(==)(varSet1::VariableSet,varSet2::VariableSet)::Bool
    if varSet1.loBs == varSet2.loBs &&
       varSet1.upBs == varSet2.upBs &&
       varSet1.vals == varSet2.vals &&
       varSet1.dscIndices == varSet2.dscIndices &&
       varSet1.bnrIndices == varSet2.bnrIndices &&
       varSet1.sos1Groups == varSet2.sos1Groups &&
       varSet1.pseudoCosts == varSet2.pseudoCosts
        return true
    end
    return false
end



function Base.:(!=)(varSet1::VariableSet,varSet2::VariableSet)::Bool
    if varSet1.loBs == varSet2.loBs &&
       varSet1.upBs == varSet2.upBs &&
       varSet1.vals == varSet2.vals &&
       varSet1.dscIndices == varSet2.dscIndices &&
       varSet1.bnrIndices == varSet2.bnrIndices &&
       varSet1.sos1Groups == varSet2.sos1Groups &&
       varSet1.pseudoCosts == varSet2.pseudoCosts
        return false
    end
    return true
end


function Base.insert!(varSet1::VariableSet,varSet2::VariableSet,insertionPoint::Int)::Nothing

    # collect info
    numNewVariables = length(varSet2.loBs)
    sos1Offset = maximum(varSet1.sos1Groups)

    numNewDiscreteVariables = length(varSet2.dscIndices)
    dscInsertionPoint = findfirst(varSet1.dscIndices.>=insertionPoint)
    if isnothing(dscInsertionPoint) dscInsertionPoint = length(varSet1.dscIndices)+1 end

    numNewBinaryVariables = length(varSet2.bnrIndices)
    bnrInsertionPoint = findfirst(varSet1.bnrIndices.>=insertionPoint)
    if isnothing(bnrInsertionPoint) bnrInsertionPoint = length(varSet1.bnrIndices)+1 end


    # append variables
    splice!(varSet1.loBs,insertionPoint:insertionPoint-1,copy(varSet2.loBs))
    splice!(varSet1.upBs,insertionPoint:insertionPoint-1,copy(varSet2.upBs))
    splice!(varSet1.vals,insertionPoint:insertionPoint-1,copy(varSet2.vals))

    splice!(varSet1.dscIndices,dscInsertionPoint:dscInsertionPoint-1,copy(varSet2.dscIndices))
    @. varSet1.dscIndices[dscInsertionPoint:dscInsertionPoint+numNewDiscreteVariables-1] += insertionPoint - 1
    @. varSet1.dscIndices[dscInsertionPoint+numNewDiscreteVariables:end] += numNewVariables

    splice!(varSet1.bnrIndices,bnrInsertionPoint:bnrInsertionPoint-1,copy(varSet2.bnrIndices))
    @. varSet1.bnrIndices[bnrInsertionPoint:bnrInsertionPoint+numNewBinaryVariables-1] += insertionPoint - 1
    @. varSet1.bnrIndices[bnrInsertionPoint+numNewBinaryVariables:end] += numNewVariables

    splice!(varSet1.sos1Groups,dscInsertionPoint:dscInsertionPoint-1,@. varSet2.sos1Groups + sos1Offset*(varSet2.sos1Groups!=0))
    varSet1.pseudoCosts = (vcat(varSet1.pseudoCosts[1][1:dscInsertionPoint-1,:],varSet2.pseudoCosts[1],varSet1.pseudoCosts[1][dscInsertionPoint:end,:]),
                                vcat(varSet1.pseudoCosts[2][1:dscInsertionPoint-1,:],varSet2.pseudoCosts[2],varSet1.pseudoCosts[2][dscInsertionPoint:end,:]))
    return

end


function Base.append!(varSet1::VariableSet,varSet2::VariableSet)::Nothing

    # collect info
    indicesOffset = get_size(varSet1)
    sos1Offset = maximum(varSet1.sos1Groups)
    # append variables
    varSet1.loBs = vcat(varSet1.loBs,varSet2.loBs)
    varSet1.upBs = vcat(varSet1.upBs,varSet2.upBs)
    varSet1.vals = vcat(varSet1.vals,varSet2.vals)
    varSet1.dscIndices = vcat(varSet1.dscIndices,@. varSet2.dscIndices + indicesOffset)
    varSet1.bnrIndices = vcat(varSet1.bnrIndices,@. varSet2.bnrIndices + indicesOffset)
    varSet1.sos1Groups = vcat(varSet1.sos1Groups,@. varSet2.sos1Groups + sos1Offset*(varSet2.sos1Groups!=0))
    varSet1.pseudoCosts = (vcat(varSet1.pseudoCosts[1],varSet2.pseudoCosts[1]),
                           vcat(varSet1.pseudoCosts[2],varSet2.pseudoCosts[2]))
    return
end

function update_bounds!(varSet::VariableSet;loBs::Vector{Float}=Float[],upBs::Vector{Float}=Float)::Nothing
    if length(loBs) > 0
        @assert length(loBs) == length(varSet.loBs) == length(varSet.upBs)
        @. varSet.loBs = loBs
    end
    if length(upBs) > 0
        @assert length(upBs) == length(varSet.loBs) == length(varSet.upBs)
        @. varSet.upBs = upBs
    end
    return
end

function update_bounds!(varSet::VariableSet,indices::Union{Vector{Int},UnitRange{Int}};loBs::Vector{Float}=Float[],upBs::Vector{Float}=Float)::Nothing
    if length(loBs) > 0
        @assert length(loBs) == length(indices)
        @. varSet.loBs[indices] = loBs
    end
    if length(upBs) > 0
        @assert length(upBs) == length(indices)
        @. varSet.upBs[indices] = upBs
    end
    return
end

function is_mixedBinary(varSet::VariableSet)::Bool
    return length(varSet.bnrIndices) == length(varSet.dscIndices)
end


## Serialization (not fundamental) used to store or send
function serialSize(varSet::VariableSet)::Int
    numVars = get_size(varSet)
    numDsc = get_numDiscrete(varSet)
    return 2 + 2*numVars + 2*numDsc
end

function serialize(varSet::VariableSet)::SerialData
    # allocate memory
    numVars = get_size(varSet)
    numDsc = get_numDiscrete(varSet)
    serial = SerialData(Vector{Float}(undef,serialSize(varSet)))

    serialize_in!(serial,varSet,offset = 0)
    return serial
end

function serialize_in!(serial::SerialData,varSet::VariableSet;offset::Int=0)::Int

    # check input
    numVars = get_size(varSet)
    numDsc = get_numDiscrete(varSet)
    @assert length(serial) >= offset + serialSize(varSet)

    # header
    serial[offset+1] = numVars
    serial[offset+2] = numDsc
    offset += 2

    # lower bounds
    serial[offset+1:offset+numVars] = varSet.loBs
    offset += numVars

    # upper bounds
    serial[offset+1:offset+numVars] = varSet.upBs
    offset += numVars

    # discrete indices
    serial[offset+1:offset+numDsc] = Vector{Float}(varSet.dscIndices)
    offset += numDsc

    # sos1 groups
    serial[offset+1:offset+numDsc] = Vector{Float}(varSet.sos1Groups)
    offset += numDsc

    return offset
end


function VariableSet(serial::SerialData;offset::Int=0)::Tuple{VariableSet,Int}

    # read header
    numVars = Int(serial[offset+1])
    numDsc = Int(serial[offset+2])
    offset += 2

    # check input
    @assert length(serial) >= offset + 2*numVars + 2*numDsc

    # build the variable set
    varSet = VariableSet(loBs=serial[offset+1:offset+numVars],
                        upBs=serial[offset+numVars+1:offset+2*numVars],
                        dscIndices=Vector{Int}(serial[offset+2*numVars+1:offset+2*numVars+numDsc]),
                        sos1Groups=Vector{Int}(serial[offset+2*numVars+numDsc+1:offset+2*numVars+2*numDsc]))

    return (varSet,offset + 2*numVars + 2*numDsc)
end


function integralize!(varSet::VariableSet,newDscIndices::Vector{Int})::Nothing
    unique!(sort!(append!(varSet.dscIndices,newDscIndices)))
    return
end
