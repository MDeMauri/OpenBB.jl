# @Author: Massimo De Mauri <massimodemauri>
# @Date:   2019-02-06T16:19:50+01:00
# @Email:  massimo.demauri@gmail.com
# @Project: OpenBB
# @Filename: problem_definitions.jl
# @Last modified by:   massimo
# @Last modified time: 2020-02-26T17:38:41+01:00
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



# serialization (Not fundamental. Used to store and send problems)
function Problem(serial::SerialData;offset::Int=0)::Tuple{Problem,Int}

    # check input
    @assert length(serial) >= offset + 2

    # header
    cnsType = Int(serial[offset+1])
    objType = Int(serial[offset+2])
    offset += 2

    (varSet,offset) = VariableSet(serial,offset=offset)

    if cnsType == 1  # LinearConstraintSet
        (cnsSet,offset) = LinearConstraintSet(serial,offset=offset)
    else
        @error "Constraint Set Type Unknown"
    end


    if objType == 1 # LinearObjective
        (objFun,offset) = LinearObjective(serial,offset=offset)
    elseif objType == 2 # QuadraticObjective
        (objFun,offset) = QuadraticObjective(serial,offset=offset)
    else
        @error "Objective Function Type Unknown"
    end

    return (Problem(varSet=varSet,cnsSet=cnsSet,objFun=objFun),offset)
end


function serial_size(problem::Problem)::Int
    return 2 + serial_size(problem.varSet) + serial_size(problem.cnsSet) + serial_size(problem.objFun)
end

function serialize(problem::Problem)::SerialData
    # assign memory
    serial = SerialData(Array{Float64,1}(undef,serial_size(problem)))
    # serialize
    serialize_in!(serial,problem,offset=0)
    return serial
end

function serialize_in!(serial::SerialData,problem::Problem;offset::Int=0)::Int

    # check input
    @assert length(serial) >= offset + 2 + serial_size(problem.varSet) + serial_size(problem.cnsSet) + serial_size(problem.objFun)

    # header
    if problem.cnsSet isa LinearConstraintSet
        serial[offset+1] = 1.0
    else
        @error "Serialization not implemented for "*str(typeof(problem.cnsSet))
    end
    if problem.objFun isa LinearObjective
        serial[offset+2] = 1.0
    elseif problem.objFun isa QuadraticObjective
        serial[offset+2] = 2.0
    else
        @error "Serialization not implemented for "*str(typeof(problem.objFun))
    end
    offset += 2

    # variable set
    offset = serialize_in!(serial,problem.varSet,offset=offset)

    # constraint set
    offset = serialize_in!(serial,problem.cnsSet,offset=offset)

    # objective function
    offset = serialize_in!(serial,problem.objFun,offset=offset)

    return offset
end
