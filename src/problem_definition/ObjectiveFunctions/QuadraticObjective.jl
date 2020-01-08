# @Author: Massimo De Mauri <massimo>
# @Date:   2019-05-23T18:34:36+02:00
# @Email:  massimo.demauri@gmail.com
# @Filename: QuadraticObjective.jl
# @Last modified by:   massimo
# @Last modified time: 2020-01-08T19:18:57+01:00
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

# copy functions
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
    return .5*point'*objective.Q*point + objective.L'*b
end

# evaluate the jacobian of the objective function in the given point
function evaluate_jacobian(objective::QuadraticObjective{T1,T2},point::Array{Float64,1})::AbstractArray where T1<:Union{Array{Float64,2},SparseMatrixCSC{Float64,Int}}  where T2<:Union{Array{Float64,1},SparseVector{Float64,Int}}
    @assert length(point) == size(objective.Q,2)
    return point'*objective.Q + objective.L'
end

# evaluate the hessian of the objective function in the given point
function evaluate_hessian(objective::QuadraticObjective{T1,T2},point::Array{Float64,1})::AbstractArray where T1<:Union{Array{Float64,2},SparseMatrixCSC{Float64,Int}}  where T2<:Union{Array{Float64,1},SparseVector{Float64,Int}}
    @assert length(point) == size(objective.Q,2)
    return objective.Q
end
