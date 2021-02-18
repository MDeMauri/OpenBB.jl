# @Author: Massimo De Mauri <massimo>
# @Date:   2020-11-06T13:59:07+01:00
# @Email:  massimo.demauri@gmail.com
# @Filename: sparsity_utils.jl
# @Last modified by:   massimo
# @Last modified time: 2021-01-14T17:30:35+01:00
# @License: LGPL-3.0
# @Copyright: {{copyright}}

export speye, spones,Infs,NaNs


# helper functions
function speye(n::Int)::SparseMatrixCSC{Float}
    indices = collect(1:n)
    return sparse(indices,indices,1.)
end

function speye(type::Type,n::Int)::SparseMatrixCSC{Float}
    indices = collect(1:n)
    return sparse(indices,indices,one(type))
end

function spones(n::Integer)::SparseVector{Float}
    return sparsevec(1:n,1.0,n)
end
function spones(type::Type,n::Integer)::SparseVector
    return sparsevec(collect(1:n),type(1),n)
end
function spones(n::Integer,m::Integer)::SparseMatrixCSC{Float}
    cols = Vector{Int}(undef,n*m)
    rows = Vector{Int}(undef,n*m)
    for k in 1:m
        @. cols[(k-1)*n+1:k*n] = 1:n
        @. rows[(k-1)*n+1:k*n] = k
    end
    return sparse(cols,rows,1.0,n,m)
end
function spones(type::Type,n::Integer,m::Integer)::SparseMatrixCSC
    cols = Vector{Int}(undef,n*m)
    rows = Vector{Int}(undef,n*m)
    for k in 1:m
        @. cols[(k-1)*n+1:k*n] = 1:n
        @. rows[(k-1)*n+1:k*n] = k
    end
    return sparse(cols,rows,type(1),n,m)
end



function findFirstNZs(A::Union{Matrix{T},SpMatrix{T}},dimension::Int)::Vector{Int} where T <: Number
    @assert dimension <= 2
    if dimension == 1
        out = Vector{Int}(undef,size(A,2))
        for k in 1:length(out)
            try
                out[k] = findfirst(!iszero,A[:,k])
            catch
                out[k] = 0
            end
        end
    else
        out = Vector{Int}(undef,size(A,1))
        for k in 1:length(out)
            try
                out[k] = findfirst(!iszero,A[k,:])
            catch
                out[k] = 0
            end
        end
    end
    return out
end


function findLastNZs(A::Union{Vector{T},SpMatrix{T}},dimension::Int)::Vector{Int} where T <: Number
    @assert dimension <= 2
    if dimension == 1
        out = Vector{Int}(undef,size(A,2))
        for k in 1:length(out)
            try
                out[k] = findlast(!iszero,A[:,k])
            catch
                out[k] = 0
            end
        end
    else
        out = Vector{Int}(undef,size(A,1))
        for k in 1:length(out)
            try
                out[k] = findlast(!iszero,A[k,:])
            catch
                out[k] = 0
            end
        end
    end
    return out
end

function sparsityMatrix(matrix::SpMatrix)::SpMatrix{Bool}
    return SpMatrix{Bool}(matrix.m,matrix.n,matrix.colptr,matrix.rowval,ones(Bool,length(matrix.nzval)))
end
function sparsityMatrix(matrix::Matrix)::SpMatrix{Bool}
    return sparsityMatrix(sparse(matrix))
end

function sparsityVector(vector::SpVector)::SpVector{Bool}
    return SpVector{Bool}(vector.n,vector.nzind,ones(Bool,length(vector.nzval)))
end
function sparsityVector(vector::Vector)::SpVector{Bool}
    return sparsityVector(sparsevec(vector))
end


function dense(matrix::SpMatrix{T})::Matrix{T} where T
    return Matrix{T}(matrix)
end
function dense(vector::SpVector{T})::Vector{T} where T
    return Vector{T}(vector)
end

function SparseArrays.sparse(template::SpMatrix{Bool},values::Vector{T})::SpMatrix{T} where T
    @assert langth(values) == nnz(template)
    return SpMatrix{T}(template.n,template.m,copy(template.colptr),copy(template.rowval),values)
end

function SparseArrays.sparse(template::SpVector{Bool},values::Vector{T})::SpVector{T} where T
    @assert langth(values) == nnz(template)
    return SpVector{T}(template.n,copy(template.nzind),values)
end

function gen_colptr(colvals::Vector{Int},numCols::Int)::Vector{Int}
   colptr = Vector{Int}(undef,numCols+1)
   numel = length(colvals)
   colptr[1] = 1
   for k in 1:numCols
       if colptr[k] <= numel
           for i in colptr[k]:numel
               if colvals[i] > k
                   colptr[k+1] = i
                   break
               elseif i==numel && colvals[i] == k
                   colptr[k+1] = i+1
               end
           end
       else
           colptr[k+1:end] .= colptr[k]
           break
       end
   end
   return colptr
end

function noncopying_sparse(rows::Vector{Int},cols::Vector{Int},vals::Vector{T},numRows::Int,numCols::Int)::SpMatrix{T} where T
    return SpMatrix{T}(numRows,numCols,gen_colptr(cols,numCols),rows,vals)
end

function noncopying_sparsevec(indices::Vector{Int},values::Vector{T},numIndices::Int)::SpVector{T} where T
    return SpVector{T}(numIndices,indices,values)
end



function spy(A::AbstractMatrix)::Nothing
    out = ""
    for i in 1:size(A,1)
       out*="|"
       for j in 1:size(A,2)
           if iszero(A[i,j])
               out *= '.'
           else
               out*='*'
           end
       end
       out*="|\n"
    end
    print(out)
    return
end
