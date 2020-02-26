# @Author: Massimo De Mauri <massimo>
# @Date:   2019-05-23T18:25:57+02:00
# @Email:  massimo.demauri@gmail.com
# @Filename: LinearConstraintSet.jl
# @Last modified by:   massimo
# @Last modified time: 2020-02-26T20:44:30+01:00
# @License: LGPL-3.0
# @Copyright: {{copyright}}


# named constructor
function LinearConstraintSet(;A::T,loBs::Array{Float64,1},upBs::Array{Float64,1})::LinearConstraintSet{T} where T<:Union{Array{Float64,2},SparseMatrixCSC{Float64,Int}}
    return LinearConstraintSet{T}(A,loBs,upBs)
end

# copy functions (Fundamental. Those are used in Branch and Bound)
import Base.copy
function copy(constraintSet::LinearConstraintSet)::LinearConstraintSet
    return LinearConstraintSet(constraintSet.A,constraintSet.loBs,constraintSet.upBs)
end

import Base.deepcopy
function deepcopy(constraintSet::LinearConstraintSet)::LinearConstraintSet
    return LinearConstraintSet(deepcopy(constraintSet.A),copy(constraintSet.loBs),copy(constraintSet.upBs))
end

# type conversion
function LinearConstraintSet{T}(constraintSet::LinearConstraintSet{T})::LinearConstraintSet{T} where T<:Union{Array{Float64,2},SparseMatrixCSC{Float64,Int}}
    return constraintSet
end

function LinearConstraintSet{T}(constraintSet::LinearConstraintSet)::LinearConstraintSet{T} where T<:Union{Array{Float64,2},SparseMatrixCSC{Float64,Int}}
    return LinearConstraintSet(T(constraintSet.A),constraintSet.loBs,constraintSet.upBs)
end

import SparseArrays.sparse
function sparse(constraintSet::LinearConstraintSet{T})::LinearConstraintSet{SparseMatrixCSC{Float64,Int}} where T<:Union{Array{Float64,2},SparseMatrixCSC{Float64,Int}}
    return LinearConstraintSet{SparseMatrixCSC{Float64,Int}}(sparse(constraintSet.A),constraintSet.loBs,constraintSet.upBs)
end


# inspect functions (Fundamental. Those are used in Branch and Bound)
function get_numVariables(constraintSet::LinearConstraintSet{T})::Int where T<:Union{Array{Float64,2},SparseMatrixCSC{Float64,Int}}
    return size(constraintSet.A,2)
end

function get_size(constraintSet::LinearConstraintSet{T})::Int where T<:Union{Array{Float64,2},SparseMatrixCSC{Float64,Int}}
    return size(constraintSet.A,1)
end

function get_bounds(constraintSet::LinearConstraintSet{T})::Tuple{Array{Float64,1},Array{Float64,1}} where T<:Union{Array{Float64,2},SparseMatrixCSC{Float64,Int}}
    return (constraintSet.loBs,constraintSet.upBs)
end


function get_sparsity(constraintSet::LinearConstraintSet{T})::Tuple{Array{Int,1},Array{Int,1}} where T<:Union{Array{Float64,2},SparseMatrixCSC{Float64,Int}}
    return findnz(sparse(constraintSet.A))[1:2]
end

function get_sparsity(constraintSet::LinearConstraintSet{T},index::Int)::Array{Int,1} where T<:Union{Array{Float64,2},SparseMatrixCSC{Float64,Int}}
    return findnz(sparse(constraintSet.A[index,:]))[1]
end

function get_firstNZs(constraintSet::LinearConstraintSet{T},dimension::Int)::Array{Int,1} where T<:Union{Array{Float64,2},SparseMatrixCSC{Float64,Int}}
    @assert 0 < dimension <= 2
    return findFirstNZs(constraintSet.A,dimension)
end

function get_firstNZs(constraintSet::LinearConstraintSet{T},indices::Array{Int,1},dimension::Int)::Array{Int,1} where T<:Union{Array{Float64,2},SparseMatrixCSC{Float64,Int}}
    @assert 0 < dimension <= 2
    if dimension == 1
        return findFirstNZs(constraintSet.A[:,indices],dimension)
    else
        return findFirstNZs(constraintSet.A[indices,:],dimension)
    end
end

function get_lastNZs(constraintSet::LinearConstraintSet{T},dimension::Int)::Array{Int,1} where T<:Union{Array{Float64,2},SparseMatrixCSC{Float64,Int}}
    @assert 0 < dimension <= 2
    return findLastNZs(constraintSet.A,dimension)
end

function get_lastNZs(constraintSet::LinearConstraintSet{T},indices::Array{Int,1},dimension::Int)::Array{Int,1} where T<:Union{Array{Float64,2},SparseMatrixCSC{Float64,Int}}
    @assert 0 < dimension <= 2
    if dimension == 1
        return findLastNZs(constraintSet.A[:,indices],dimension)
    else
        return findLastNZs(constraintSet.A[indices,:],dimension)
    end
end


function get_linearConstraints(constraintSet::LinearConstraintSet{T})::T where T<:Union{Array{Float64,2},SparseMatrixCSC{Float64,Int}}
    return constraintSet.A
end

# update functions (Not fundamental. Those are used only in updating the problem)
function update_bounds!(constraintSet::LinearConstraintSet;loBs::Array{Float64,1}=Float64[],upBs::Array{Float64,1}=Float64[])::Nothing
    if length(loBs) > 0
        @assert length(loBs) == length(constraintSet.loBs) == length(constraintSet.upBs)
        @. constraintSet.loBs = loBs
    end
    if length(upBs) > 0
        @assert length(upBs) == length(constraintSet.loBs) == length(constraintSet.upBs)
        @. constraintSet.upBs = upBs
    end
    return
end

function update_bounds!(constraintSet::LinearConstraintSet,indices::Array{Int,1};loBs::Array{Float64,1}=Float64[],upBs::Array{Float64,1}=Float64[])::Nothing
    if length(loBs) > 0
        @assert length(loBs) == length(indices)
        @. constraintSet.loBs[indices] = loBs
    end
    if length(upBs) > 0
        @assert length(upBs) == length(indices)
        @. constraintSet.upBs[indices] = upBs
    end
    return
end

function remove_variables!(constraintSet::LinearConstraintSet{T},indices::Array{Int,1})::Nothing where T<:Union{Array{Float64,2},SparseMatrixCSC{Float64,Int}}
    toKeep = filter(x->!(x in indices), collect(1:get_numVariables(constraintSet)))
    constraintSet.A = constraintSet.A[:,toKeep]
    return
end

function insert_variables!(constraintSet::LinearConstraintSet{T},numVariables::Int,insertionPoint::Int)::Nothing where T<:Union{Array{Float64,2},SparseMatrixCSC{Float64,Int}}
    @assert numVariables>=0
    @assert 0<insertionPoint<=get_numVariables(constraintSet)+1
    constraintSet.A = hcat(constraintSet.A[:,1:insertionPoint-1],zeros(size(constraintSet.A,1),numVariables),constraintSet.A[:,insertionPoint:end])
    return
end

function append_variables!(constraintSet::LinearConstraintSet{T},numVariables::Int)::Nothing where T<:Union{Array{Float64,2},SparseMatrixCSC{Float64,Int}}
    insert_variables!(constraintSet,numVariables,get_numVariables(constraintSet)+1)
    return
end

function remove_constraints!(constraintSet::LinearConstraintSet{T},indices::Array{Int,1})::Nothing  where T<:Union{Array{Float64,2},SparseMatrixCSC{Float64,Int}}
    toKeep = filter(x->!(x in indices), collect(1:get_size(constraintSet)))
    constraintSet.A = constraintSet.A[toKeep,:]
    constraintSet.loBs = constraintSet.loBs[toKeep]
    constraintSet.upBs = constraintSet.upBs[toKeep]
    return
end

import Base.insert!
function insert!(constraintSet1::LinearConstraintSet{T1},constraintSet2::LinearConstraintSet{T2},insertionPoint::Int)::Nothing where T1<:Union{Array{Float64,2},SparseMatrixCSC{Float64,Int}} where T2<:Union{Array{Float64,2},SparseMatrixCSC{Float64,Int}}
    constraintSet1.A = vcat(constraintSet1.A[1:insertionPoint-1,:],constraintSet2.A,constraintSet1.A[insertionPoint:end,:])
    splice!(constraintSet1.loBs,insertionPoint:insertionPoint-1,copy(constraintSet2.loBs))
    splice!(constraintSet1.upBs,insertionPoint:insertionPoint-1,copy(constraintSet2.upBs))
    return
end


import Base.append!
function append!(constraintSet1::LinearConstraintSet{T1},constraintSet2::LinearConstraintSet{T2})::Nothing where T1<:Union{Array{Float64,2},SparseMatrixCSC{Float64,Int}} where T2<:Union{Array{Float64,2},SparseMatrixCSC{Float64,Int}}
    return insert!(constraintSet1,constraintSet2,get_size(constraintSet1)+1)
end

function concat(constraintSet1::LinearConstraintSet{T1},constraintSet2::LinearConstraintSet{T2})::LinearConstraintSet where T1<:Union{Array{Float64,2},SparseMatrixCSC{Float64,Int}} where T2<:Union{Array{Float64,2},SparseMatrixCSC{Float64,Int}}

    return LinearConstraintSet(vcat(constraintSet1.A,constraintSet2.A),
                               vcat(constraintSet1.loBs,constraintSet2.loBs),
                               vcat(constraintSet1.upBs,constraintSet2.upBs))
end

function permute_constraints!(constraintSet::LinearConstraintSet{T},permutation::Array{Int,1})::Nothing where T<:Union{Array{Float64,2},SparseMatrixCSC{Float64,Int}}
    constraintSet.A = constraintSet.A[permutation,:]
    constraintSet.loBs = constraintSet.loBs[permutation]
    constraintSet.upBs = constraintSet.upBs[permutation]
    return
end

# give the value of the constraints in the given point
function evaluate(constraintSet::LinearConstraintSet{T},point::Array{Float64,1})::Array{Float64,1} where T<:Union{Array{Float64,2},SparseMatrixCSC{Float64,Int}}
    @assert length(point) == size(constraintSet.A,2)
    return constraintSet.A*point
end

# evaluate the jacobian of the constraints in the given point
function evaluate_jacobian(constraintSet::LinearConstraintSet{T},point::Array{Float64,1})::AbstractMatrix  where T<:Union{Array{Float64,1},SparseVector{Float64,Int}}
    @assert length(point) == length(objective.L)
    return copy(constraintSet.A)
end


# Serialization (not fundamental) used to store or to send
function LinearConstraintSet(serial::SerialData;offset::Int=0)::Tuple{LinearConstraintSet,Int}

    # header
    numCnss = Int(serial[offset+1])
    numVars = Int(serial[offset+2])
    numNZs = Int(serial[offset+3])
    offset += 3


    if numNZs == -1 # dense constraint matrix
        # check input
        @assert length(serial) >= offset + numVars*numCnss + 2*numCnss

        # reconstruct constraint matrix
        A = zeros(numCnss,numVars)
        for k in 1:numCnss
            A[k,:] = serial[offset+1:offset+numVars]
            offset += numVars
        end
    else # sparse constraint matrix

        # check input
        @assert length(serial) >= offset + 3*numNZs + 2*numCnss

        # reconstruct constraint matrix
        A = sparse(Array{Int,1}(serial[offset+1:offset+numNZs]),
                   Array{Int,1}(serial[offset+numNZs+1:offset+2*numNZs]),
                   serial[offset+2*numNZs+1:offset+3*numNZs],
                   numCnss,numVars)
        offset += 3*numNZs
    end

    # bounds
    loBs = serial[offset+1:offset+numCnss]
    offset += numCnss
    upBs = serial[offset+1:offset+numCnss]
    offset += numCnss


    # reconstruct the constraint set
    return (LinearConstraintSet(A,loBs,upBs),offset)

end


function serial_size(cnsSet::LinearConstraintSet{T})::Int where T<:Union{Array{Float64,2},SparseMatrixCSC{Float64,Int}}
    # estimate necessary memory
    numCnss = size(cnsSet.A,1)
    numVars = size(cnsSet.A,2)
    if T == Array{Float64,2}
        return 3 + numCnss*numVars + 2*numCnss
    else
        return 3 + 3*nnz(cnsSet.A) + 2*numCnss
    end
end

function serialize(cnsSet::LinearConstraintSet{T})::SerialData where T<:Union{Array{Float64,2},SparseMatrixCSC{Float64,Int}}

    # allocate memory
    serial = SerialData(Array{Float64,1}(undef,serial_size(cnsSet)))

    # write serialized data
    serialize_in!(serial,cnsSet,offset=0)
    return serial
end


function serialize_in!(serial::SerialData,cnsSet::LinearConstraintSet{T};offset::Int=0)::Int where T<:Union{Array{Float64,2},SparseMatrixCSC{Float64,Int}}

    # check input
    numCnss = size(cnsSet.A,1)
    numVars = size(cnsSet.A,2)

    if T == Array{Float64,2}
        # check size
        @assert length(serial) >= offset + 3 + numCnss*numVars + 2*numCnss

        # header
        serial[offset+1] = numCnss
        serial[offset+2] = numVars
        serial[offset+3] = -1
        offset += 3

        # constraint matrix
        for k in 1:numCnss
            serial[offset+1:offset+numVars] = cnsSet.A[k,:]
            offset += numVars
        end
    else
        # check size
        numNZs = nnz(cnsSet.A)
        @assert length(serial) >= offset + 3 + 3*numNZs + 2*numCnss

        # header
        serial[offset+1] = numCnss
        serial[offset+2] = numVars
        serial[offset+3] = numNZs
        offset += 3

        # constraint matrix
        data_ = findnz(cnsSet.A)
        for k in 1:3
            serial[offset+(k-1)*numNZs+1:offset+k*numNZs] = data_[k]
        end
        offset += 3*numNZs
    end

    # bounds
    serial[offset+1:offset+numCnss] = cnsSet.loBs
    offset += numCnss
    serial[offset+1:offset+numCnss] = cnsSet.upBs
    offset += numCnss

    return offset
end
