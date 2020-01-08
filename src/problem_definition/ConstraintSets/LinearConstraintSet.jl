# @Author: Massimo De Mauri <massimo>
# @Date:   2019-05-23T18:25:57+02:00
# @Email:  massimo.demauri@gmail.com
# @Filename: LinearConstraintSet.jl
# @Last modified by:   massimo
# @Last modified time: 2020-01-08T19:21:45+01:00
# @License: LGPL-3.0
# @Copyright: {{copyright}}

# constructirs and copy functions (Fundamental. Those are used in Branch and Bound)
# named constructor
function LinearConstraintSet(;A::T,loBs::Array{Float64,1},upBs::Array{Float64,1})::LinearConstraintSet{T} where T<:Union{Array{Float64,2},SparseMatrixCSC{Float64,Int}}
    return LinearConstraintSet{T}(A,loBs,upBs)
end


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

# returns the number of nonzeros in the constraints matrix
import SparseArrays.nnz
function nnz(constraintSet::LinearConstraintSet{T})::Int where T<:Union{Array{Float64,2},SparseMatrixCSC{Float64,Int}}
    return nnz(sparse(constraintSet.A))
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
