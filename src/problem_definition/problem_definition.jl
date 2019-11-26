# @Author: Massimo De Mauri <massimodemauri>
# @Date:   2019-02-06T16:19:50+01:00
# @Email:  massimo.demauri@gmail.com
# @Project: OpenBB
# @Filename: problem_definitions.jl
# @Last modified by:   massimo
# @Last modified time: 2019-11-22T12:49:02+01:00
# @License: LGPL-3.0
# @Copyright: {{copyright}}


# load the possible components definitions
include("./ObjectiveFunctions/ObjectiveFunctions.jl")
include("./ConstraintSets/ConstraintSets.jl")
include("./VariableSet.jl")


# abstract and null types
abstract type AbstractProblem end
struct NullProblem <: AbstractProblem end

# concrete problem
mutable struct Problem{T1<:AbstractObjective,T2<:AbstractConstraintSet}
    objFun::T1
    cnsSet::T2
    varSet::VariableSet
end

# constructors and copy functions (Fundamental. These are used by branch and bound.)
function Problem(;objFun::AbstractObjective,cnsSet::AbstractConstraintSet,varSet::VariableSet)::Problem
    @assert get_numVariables(objFun) == get_numVariables(cnsSet) == get_size(varSet)
    return Problem(objFun,cnsSet,varSet)
end

import Base.copy
function copy(problem::Problem)::Problem
    return Problem(copy(problem.objFun),copy(problem.cnsSet),copy(problem.varSet))
end

import Base.deepcopy
function deepcopy(problem::Problem)::Problem
    return Problem(deepcopy(problem.objFun),deepcopy(problem.cnsSet),deepcopy(problem.varSet))
end

import SparseArrays.sparse
function sparse(problem::Problem)::Problem
    return Problem(sparse(problem.objFun),sparse(problem.cnsSet),problem.varSet)
end


# inspect functions (Fundamental. These are used by branch and bound.)
function get_numVariables(problem::Problem)::Int
    return get_size(problem.varSet)
end

function get_numDiscreteVariables(problem::Problem)::Int
    return get_numDiscrete(problem.varSet)
end

function get_variableBounds(problem::Problem)::Tuple{Array{Float64,1},Array{Float64,1}}
    return get_bounds(problem.varSet)
end

function get_discreteIndices(problem::Problem)::Array{Int,1}
    return get_discreteIndices(problem.varSet)
end

function get_sos1Groups(problem)::Array{Int,1}
    return get_sos1Groups(problem.varSet)
end

function get_pseudoCosts(problem::Problem)::Tuple{Array{Float64,2},Array{Int,2}}
    return get_pseudoCosts(problem.varSet)
end

function get_numConstraints(problem::Problem)::Int
    return get_size(problem.cnsSet)
end

function get_constraintBounds(problem::Problem)::Tuple{Array{Float64,1},Array{Float64,1}}
    return get_bounds(problem.cnsSet)
end

function get_objective_sparsity(problem::Problem)::Any
    return get_sparsity(problem.objFun)
end

function get_constraints_sparsity(problem::Problem)::Any
    return get_sparsity(problem.cnsSet)
end

function get_constraint_sparsity(problem::Problem,index::Int)::Any
    return get_sparsity(problem.cnsSet,index)
end
