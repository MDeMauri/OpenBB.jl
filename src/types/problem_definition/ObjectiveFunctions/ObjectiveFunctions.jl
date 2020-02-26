# @Author: Massimo De Mauri <massimo>
# @Date:   2019-02-26T15:03:06+01:00
# @Email:  massimo.demauri@gmail.com
# @Filename: ObjectiveFunctions.jl
# @Last modified by:   massimo
# @Last modified time: 2019-08-27T18:14:17+02:00
# @License: LGPL-3.0
# @Copyright: {{copyright}}


# abstract and null types
abstract type AbstractObjective end
struct NullObjective <: AbstractObjective end


# linear objective
mutable struct LinearObjective{T} <: AbstractObjective where T<:Union{Array{Float64,1},SparseVector{Float64,Int}}
    L::T
end

# quadratic objective
mutable struct QuadraticObjective{T1,T2} <: AbstractObjective where T1<:Union{Array{Float64,2},SparseMatrixCSC{Float64,Int}} where T2<:Union{Array{Float64,1},SparseVector{Float64,Int}}
    Q::T1
    L::T2
end


include("./LinearObjective.jl")
include("./QuadraticObjective.jl")
