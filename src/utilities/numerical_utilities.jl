# @Author: Massimo De Mauri <massimo>
# @Date:   2019-03-01T20:22:38+01:00
# @Email:  massimo.demauri@gmail.com
# @Filename: numerical_utilities.jl
# @Last modified by:   massimo
# @Last modified time: 2019-11-20T15:29:41+01:00
# @License: LGPL-3.0
# @Copyright: {{copyright}}

const IntInf = 2^63-1



export speye, spones,Infs,NaNs

# helper functions
function speye(n::Int)::SparseMatrixCSC
    indices = collect(1:n)
    return sparse(indices,indices,1.)
end

function spones(n::Integer)::SparseMatrixCSC
    return sparse(ones(n))
end
function spones(n::Integer,m::Integer)::SparseMatrixCSC
    return sparse(ones(n,m))
end

function Infs(n::Integer)::Array{Float64,1}
    return Inf*ones(n)
end

function Infs(n::Integer,m::Integer)::Array{Float64,2}
    return Inf*ones(n,m)
end

function NaNs(n::Integer)::Array{Float64,1}
    return NaN*ones(n)
end

function NaNs(n::Integer,m::Integer)::Array{Float64,2}
    return NaN*ones(n,m)
end


function threshold(number::T1,threshold::T2)::T1 where T1 <: Real where T2 <: Real
    if number >= threshold
        return number
    else
        return T2(0)
    end
end

function clip(number::T1,threshold::T2)::T1 where T1 <: Real where T2 <: Real
    if number <= threshold
        return number
    else
        return threshold
    end
end



function findFirstNZs(A::Union{Array{T,1},SparseMatrixCSC{T,Int}},dimension::Int)::Array{Int,1} where T <: Number
    @assert dimension <= 2
    if dimension == 1
        out = Array{Int,1}(undef,size(A,2))
        for k in 1:length(out)
            try
                out[k] = findfirst(!iszero,A[:,k])
            catch
                out[k] = -1
            end
        end
    else
        out = Array{Int,1}(undef,size(A,1))
        for k in 1:length(out)
            try
                out[k] = findfirst(!iszero,A[k,:])
            catch
                out[k] = -1
            end
        end
    end
    return out
end


function findLastNZs(A::Union{Array{T,1},SparseMatrixCSC{T,Int}},dimension::Int)::Array{Int,1} where T <: Number
    @assert dimension <= 2
    if dimension == 1
        out = Array{Int,1}(undef,size(A,2))
        for k in 1:length(out)
            try
                out[k] = findlast(!iszero,A[:,k])
            catch
                out[k] = -1
            end
        end
    else
        out = Array{Int,1}(undef,size(A,1))
        for k in 1:length(out)
            try
                out[k] = findlast(!iszero,A[k,:])
            catch
                out[k] = -1
            end
        end
    end
    return out
end
