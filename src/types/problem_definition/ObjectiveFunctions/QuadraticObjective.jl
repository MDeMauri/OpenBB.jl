# @Author: Massimo De Mauri <massimo>
# @Date:   2019-05-23T18:34:36+02:00
# @Email:  massimo.demauri@gmail.com
# @Filename: QuadraticObjective.jl
# @Last modified by:   massimo
# @Last modified time: 2020-02-26T21:21:54+01:00
# @License: LGPL-3.0
# @Copyright: {{copyright}}


# constructors and copy functions (Fundamental. These are used in Branch and Bound)
# named constructor
function QuadraticObjective(;Q::T1,L::T2)::QuadraticObjective{T1,T2} where T1<:Union{Array{Float64,2},SparseMatrixCSC{Float64,Int}}  where T2<:Union{Array{Float64,1},SparseVector{Float64,Int}}
    return QuadraticObjective(Q,L)
end

# type conversions
function QuadraticObjective{T1,T2}(objective::LinearObjective{T2})::QuadraticObjective{T1,T2} where T1<:Union{Array{Float64,2},SparseMatrixCSC{Float64,Int}}  where T2<:Union{Array{Float64,1},SparseVector{Float64,Int}}
    return QuadraticObjective(T1(zeros(length(objective.L),length(objective.L))),T2(objective.L))
end

function QuadraticObjective{T1,T2}(objective::QuadraticObjective{T1,T2})::QuadraticObjective{T1,T2} where T1<:Union{Array{Float64,2},SparseMatrixCSC{Float64,Int}}  where T2<:Union{Array{Float64,1},SparseVector{Float64,Int}}
    return objective
end

function QuadraticObjective{T1,T2}(objective::QuadraticObjective)::QuadraticObjective{T1,T2} where T1<:Union{Array{Float64,2},SparseMatrixCSC{Float64,Int}}  where T2<:Union{Array{Float64,1},SparseVector{Float64,Int}}
    return QuadraticObjective(T1(objective.Q),T2(objective.L))
end

# copy functions (Fundamental)
import Base.copy
function copy(objective::QuadraticObjective)::QuadraticObjective
    return QuadraticObjective(objective.Q,objective.L)
end
import Base.deepcopy
function deepcopy(objective::QuadraticObjective)::QuadraticObjective
    return QuadraticObjective(deepcopy(objective.Q),deepcopy(objective.L))
end

import SparseArrays.sparse
function sparse(objective::QuadraticObjective)::QuadraticObjective
    return QuadraticObjective(sparse(objective.Q),sparse(objective.L))
end


# inspect functions  (Fundamental. These are used in Branch and Bound)
function get_numVariables(objective::QuadraticObjective)::Int
    return size(objective.Q,1)
end

function get_sparsity(objective::QuadraticObjective)::Tuple{Tuple{Array{Int,1},Array{Int,1}},Array{Int,1}}
    return (findnz(sparse(objective.Q))[1:2],findnz(sparse(objective.L))[1])
end

# update functions (Not fundamental These are used only for problem update)
function insert_variables!(objective::QuadraticObjective,numVariables::Int,insertionPoint::Int)::Nothing
    @assert numVariables >= 0
    @assert 0<=insertionPoint<=get_numVariables(objective)+1
    objective.Q = vcat(hcat(objective.Q[1:insertionPoint-1,1:insertionPoint-1],zeros(insertionPoint-1,numVariables),objective.Q[1:insertionPoint-1,insertionPoint:end]),
                       zeros(numVariables,numVariables+get_numVariables(objective)),
                       hcat(objective.Q[insertionPoint:end,1:insertionPoint-1],zeros(get_numVariables(objective)-insertionPoint+1,numVariables),objective.Q[insertionPoint:end,insertionPoint:end]))
    splice!(objective.L,insertionPoint:insertionPoint-1,zeros(numVariables,1))
    return
end

function append_variables!(objective::QuadraticObjective,numVariables::Int)::Nothing
    insert_variables!(objective,numVariables,get_numVariables(objective)+1)
    return
end


function remove_variables!(objective::QuadraticObjective,indices::Array{Int,1})::Nothing
    toKeep = filter(x->!(x in indices), collect(1:get_numVariables(objective)))
    objective.Q = objective.Q[toKeep,toKeep]
    objective.L = objective.L[toKeep]
    return
end


import Base.+
function +(objective1::QuadraticObjective,objective2::QuadraticObjective)::QuadraticObjective where T<:AbstractObjective
    @assert get_numVariables(objective1) == get_numVariables(objective2)
    return QuadraticObjective(objective1.Q+objective2.Q,objective1.L+objective2.L)
end

function +(objective1::QuadraticObjective,objective2::T)::QuadraticObjective where T<:AbstractObjective

    if objective2 isa NullObjective
        return objective1
    else
        return objective1 + QuadraticObjective(objective2)
    end
end



function add!(objective1::QuadraticObjective,objective2::QuadraticObjective)::Nothing
    @assert get_numVariables(objective1) == get_numVariables(objective2)
    @. objective1.Q += objective2.Q
    @. objective1.L += objective2.L
    return
end

function add!(objective1::QuadraticObjective,objective2::T)::Nothing where T<:AbstractObjective
    if objective2 isa NullObjective
        return
    else
        add!(objective1,QuadraticObjective(objective2))
        return
    end
end



# evaluate the objective function in the given point
function evaluate(objective::QuadraticObjective{T1,T2},point::Array{Float64,1})::Float64 where T1<:Union{Array{Float64,2},SparseMatrixCSC{Float64,Int}}  where T2<:Union{Array{Float64,1},SparseVector{Float64,Int}}
    @assert length(point) == size(objective.Q,2)
    return .5*point'*objective.Q*point + objective.L'*point
end

# evaluate the jacobian of the objective function in the given point
function evaluate_jacobian(objective::QuadraticObjective{T1,T2},point::Array{Float64,1})::AbstractArray where T1<:Union{Array{Float64,2},SparseMatrixCSC{Float64,Int}}  where T2<:Union{Array{Float64,1},SparseVector{Float64,Int}}
    @assert length(point) == size(objective.Q,2)
    return point'*objective.Q + objective.L'
end

# evaluate the hessian of the objective function in the given point
function evaluate_hessian(objective::QuadraticObjective{T1,T2},point::Array{Float64,1})::AbstractArray where T1<:Union{Array{Float64,2},SparseMatrixCSC{Float64,Int}}  where T2<:Union{Array{Float64,1},SparseVector{Float64,Int}}
    @assert length(point) == size(objective.Q,2)
    return copy(objective.Q)
end



# Serialization (Not fundamental. Used to store and to send)
function QuadraticObjective(serial::SerialData;offset::Int=0)::Tuple{QuadraticObjective,Int}
    # header
    numVars = Int(serial[offset+1])
    numNZsQ = Int(serial[offset+2])
    numNZsL = Int(serial[offset+3])
    offset += 3

    if numNZsQ == -1 # dense quadratic term
        # check input
        @assert length(serial) >= offset + numVars^2
        # reconstruct constraint matrix
        Q = zeros(numVars,numVars)
        for k in 1:numVars
            Q[k,:] = serial[offset+1:offset+numVars]
            offset += numVars
        end
    else # sparse constraint matrix
        # check input
        @assert length(serial) >= offset + 3*numNZsQ
        # reconstruct constraint matrix
        Q = sparse(Array{Int,1}(serial[offset+1:offset+numNZsQ]),
                   Array{Int,1}(serial[offset+numNZsQ+1:offset+2*numNZsQ]),
                   serial[offset+2*numNZsQ+1:offset+3*numNZsQ],
                   numVars,numVars)
        offset += 3*numNZsQ
    end

    # linear term
    if numNZsL == -1 # dense linear term
        @assert length(serial) >= offset + numVars
        L = serial[offset+1:offset+numVars]
        offset+=numVars
    else
        L = sparsevec(Array{Int,1}(serial[offset+1:offset+numNZsL]),
                      serial[offset+numNZsL+1:offset+2*numNZsL],
                      numVars)
        offset += 2*numNZsL
    end
    
    # reconstruct the objective function
    return (QuadraticObjective(Q,L),offset)

end

function serial_size(objFun::QuadraticObjective{T1,T2})::Int where T1<:Union{Array{Float64,2},SparseMatrixCSC{Float64,Int}} where T2<:Union{Array{Float64,1},SparseVector{Float64,Int}}
    # estimate necessary memory
    numVars = size(objFun.Q,1)
    sSize = 3
    if T1 == Array{Float64,2}
        sSize += numVars^2
    else
        sSize += 3*nnz(objFun.Q)
    end

    if T2 == Array{Float64,1}
        sSize += numVars
    else
        sSize += 2*nnz(objFun.L)
    end

    return sSize
end

function serialize(objFun::QuadraticObjective{T1,T2})::SerialData where T1<:Union{Array{Float64,2},SparseMatrixCSC{Float64,Int}} where T2<:Union{Array{Float64,1},SparseVector{Float64,Int}}

    # allocate memory
    serial = SerialData(Array{Float64,1}(undef,serial_size(objFun)))

    # write serialized data
    serialize_in!(serial,objFun,offset=0)
    return serial
end

function serialize_in!(serial::SerialData,objFun::QuadraticObjective{T1,T2};offset::Int=0)::Int where T1<:Union{Array{Float64,2},SparseMatrixCSC{Float64,Int}} where T2<:Union{Array{Float64,1},SparseVector{Float64,Int}}
    # check input
    @assert length(serial) >= offset + serial_size(objFun)

    # collect info
    numVars = size(objFun.Q,1)

    # header
    serial[offset+1] = numVars
    if T1 == Array{Float64,2}
        serial[offset+2] = -1.0
    else
        serial[offset+2] = nnz(objFun.Q)
    end

    if T2 == Array{Float64,1}
        serial[offset+3] = -1.0
    else
        serial[offset+3] = nnz(objFun.L)
    end
    offset += 3

    if T1 == Array{Float64,2} # dense quadratic term
        for k in 1:numVars
            serial[offset+1:offset+numVars] = objFun.Q[k,:]
            offset += numVars
        end
    else # sparse quadratic term
        data_ = findnz(objFun.Q)
        numNZs_ = length(data_[1])
        for k in 1:3
            serial[offset+(k-1)*numNZs_+1:offset+k*numNZs_] = data_[k]
        end
        offset += 3*numNZs_
    end

    if T2 == Array{Float64,1} #dense linear term
        serial[offset+1:offset+numVars] = objFun.L
        offset += numVars
    else # sparse linear term
        data_ = findnz(objFun.L)
        numNZs_ = length(data_[1])
        for k in 1:2
            serial[offset+(k-1)*numNZs_+1:offset+k*numNZs_] = data_[k]
        end
        offset += 2*numNZs_
    end

    return offset
end
