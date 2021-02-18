# @Author: Massimo De Mauri <massimo>
# @Date:   2019-05-23T18:33:51+02:00
# @Email:  massimo.demauri@gmail.com
# @Filename: LinearObjective.jl
# @Last modified by:   massimo
# @Last modified time: 2021-01-19T23:36:30+01:00
# @License: LGPL-3.0
# @Copyright: {{copyright}}

function is_linear(objFun::LinearObjective)::Bool
    return true
end
function is_quadratic(objFun::LinearObjective)::Bool
    return true
end
function is_convex(objFun::LinearObjective)::Bool
    return true
end


# constructors and copy functions (Fundamental. These are used in Branch and Bound)
# named constructor
function LinearObjective(;L::T,c::Float=0.0)::LinearObjective{T} where T<:Union{Vector{Float},SpVector{Float}}
    return LinearObjective(L,c)
end

function Base.copy(objFun::LinearObjective)::LinearObjective
    return LinearObjective(objFun.L,objFun.c)
end

function Base.deepcopy(objFun::LinearObjective)::LinearObjective
    return LinearObjective(deepcopy(objFun.L),objFun.c)
end

# type conversions
function LinearObjective{T}(objFun::LinearObjective{T})::LinearObjective{T} where T<:Union{Vector{Float},SpVector{Float}}
    return deepcopy(objFun)
end

function LinearObjective{T}(objFun::LinearObjective)::LinearObjective{T} where T<:Union{Vector{Float},SpVector{Float}}
    return LinearObjective(T(objFun.L),objFun.c)
end

function LinearObjective{T}(objFun::QuadraticObjective)::LinearObjective{T} where T<:Union{Vector{Float},SpVector{Float}}
    @assert iszero(objFun.Q)
    return LinearObjective(T(objFun.L),objFun.c)
end

function SparseArrays.sparse(objFun::LinearObjective)::LinearObjective{SpVector{Float}}
    return LinearObjective{SpVector{Float}}(objFun)
end

# inspect functions (Fundamental. These are used in Branch and Bound)
function get_numVariables(objFun::LinearObjective)::Int
    return length(objFun.L)
end

function get_dependency(objFun::LinearObjective)::Vector{Int}
    if objFun.L isa Vector
        return findall(!iszero,objFun.L)
    else
        return findnz(objFun.L)[1]
    end
end

function get_gradientSparsity(objFun::LinearObjective)::Vector{Int}
    if objFun.L isa Vector
        return findall(!iszero,objFun.L)[1]
    else
        return findnz(objFun.L)[1]
    end
end

function get_hessianSparsity(objFun::LinearObjective)::Tuple{Vector{Int},Vector{Int}}
    return (Int[],Int[])
end

function get_gradientType(objFun::LinearObjective{T})::Type where T<:Union{Vector{Float},SpVector{Float}}
    return T
end

function get_hessianType(objFun::LinearObjective{T})::Type where T<:Union{Vector{Float},SpVector{Float}}
    if T <: Vector
        return Matrix{Float}
    else
        return SpMatrix{Float}
    end
end

# evaluate the objFun function in the given point
function evaluate(objFun::LinearObjective,point::VirtualVector{Float})::Float
    @assert length(point) == length(objFun.L)
    return objFun.L'*point + objFun.c
end

# evaluate the gradient of the objFun function in the given point
function evaluate_gradient(objFun::LinearObjective{TL},point::VirtualVector{Float})::TL where TL<:Union{Vector{Float},SpVector{Float}}
    @assert length(point) == length(objFun.L)
    return copy(objFun.L)
end

# evaluate the hessian of the objFun function in the given point
function evaluate_hessian(objFun::LinearObjective{Vector{Float}},point::VirtualVector{Float})::Matrix{Float}
    @assert length(point) == length(objFun.L)
    return zeros(length(point),length(point))
end

function evaluate_hessian(objFun::LinearObjective{SpVector{Float}},point::VirtualVector{Float})::SpMatrix{Float}
    @assert length(point) == length(objFun.L)
    return spzeros(length(point),length(point))
end


# update functions (Not Fundamental. These are used only during problem update)
function insert_variables!(objFun::LinearObjective,numVariables::Int,insertionPoint::Int)::Nothing
    @assert numVariables>=0
    @assert 0<=insertionPoint<=get_numVariables(objFun)+1
    if objFun isa LinearObjective{Vector{Float}}
        splice!(objFun.L,insertionPoint:insertionPoint-1,zeros(numVariables))
    else
        objFun.L = vcat(objFun.L[1:insertionPoint-1],spzeros(numVariables),objFun.L[insertionPoint:end])
    end
    return
end

function append_variables!(objFun::LinearObjective,numVariables::Int)::Nothing
    insert_variables!(objFun,numVariables,get_numVariables(objFun)+1)
    return
end


function remove_variables!(objFun::LinearObjective,indices::Union{Vector{Int},UnitRange{Int}})::Nothing
    toKeep = filter(x->!(x in indices), collect(1:get_numVariables(objFun)))
    objFun.L = objFun.L[toKeep]
    return
end

function fix_variables!(objFun::LinearObjective,indices::Union{Vector{Int},UnitRange{Int}},values::Vector{Float};removeFixedVariables::Bool=false)::Nothing
    toKeep = filter(x->!(x in indices), collect(1:get_numVariables(objFun)))
    objFun.c = objFun.L[indices]'*values
    if removeFixedVariables
        objFun.L = objFun.L[toKeep]
    else
        if objFun.L isa Vector
            objFun.L[indices] = zeros(length(indices))
        else
            objFun.L[indices] = spzeros(length(indices))
        end
    end
    return
end


import Base.+
function +(objFun1::LinearObjective{T},objFun2::LinearObjective{T})::LinearObjective where T<:Union{Vector{Float},SpVector{Float}}
    @assert get_numVariables(objFun1) == get_numVariables(objFun2)
    return LinearObjective(objFun1.L+objFun2.L,objFun1.c+objFun2.c)
end

function add!(objFun1::LinearObjective,objFun2::LinearObjective)::Nothing
    @assert get_numVariables(objFun1) == get_numVariables(objFun2)
    @. objFun1.L += objFun2.L
    objFun1.c += objFun2.c
    return
end


# Serialization (not fundamental) used to store or to send
function LinearObjective(serial::SerialData;offset::Int=0)::Tuple{LinearObjective,Int}
    # header
    numVars = Int(serial[offset+1])
    numNZs = Int(serial[offset+2])
    offset += 2
    # data
    if numNZs == -1 # Dense L
        @assert length(serial) >= offset + numVars + 1
        L = Vector{Float}(serial[offset+1:offset+numVars]); offset += numVars
        c = serial[offset+1]; offset += 1
        return  (LinearObjective(L,c),offset)
    else # sparse L
        @assert length(serial) >= offset + 2*numNZs + 1
        inds = Vector{Int}(serial[offset+1:offset+numNZs]); offset += numNZs
        vals = serial[offset+1:offset+numNZs]; offset += numNZs
        c = serial[offset+1]; offset += 1
        return (LinearObjective(sparsevec(inds,vals,numVars),c),offset)
    end
end


function serialSize(objFun::LinearObjective)::Int
    # estimate necessary memory
    if objFun isa LinearObjective{Vector{Float}}
        return length(objFun.L) + 3
    else
        return 2*nnz(objFun.L) + 3
    end
end


function serialize(objFun::LinearObjective)::SerialData

    # allocate memory
    serial = SerialData(Vector{Float}(undef,serialSize(objFun)))

    # write serialized data
    serialize_in!(serial,objFun,offset=0)
    return serial
end


function serialize_in!(serial::SerialData,objFun::LinearObjective;offset::Int=0)::Int

    # check input
    @assert length(serial) >= serialSize(objFun) + offset

    # store info
    numVars = get_numVariables(objFun)
    if objFun isa LinearObjective{Vector{Float}}
        numNZs = -1

        # header
        serial[offset+1] = numVars
        serial[offset+2] = numNZs
        offset += 2

        # data
        serial[offset+1:offset+numVars] = objFun.L; offset += numVars
        serial[offset+1] = objFun.c; offset += 1
    else
        numNZs = nnz(objFun.L)

        # header
        serial[offset+1] = numVars
        serial[offset+2] = numNZs
        offset += 2

        # data
        (inds,vals) = findnz(objFun.L)
        serial[offset+1:offset+numNZs] = inds; offset += numNZs
        serial[offset+1:offset+numNZs] = vals; offset += numNZs
        serial[offset+1] = objFun.c; offset += 1
    end

    return offset
end


for type in LinearObjectiveLeafTypes @assert precompile(is_linear,(type,)) end
for type in LinearObjectiveLeafTypes @assert precompile(is_quadratic,(type,)) end
for type in LinearObjectiveLeafTypes @assert precompile(is_convex,(type,)) end
for type in LinearObjectiveLeafTypes @assert precompile(copy,(type,)) end
for type in LinearObjectiveLeafTypes @assert precompile(deepcopy,(type,)) end
for type1 in LinearObjectiveLeafTypes
    for type2 in LinearObjectiveLeafTypes @assert precompile(type1,(type2,)) end
    for type2 in QuadraticObjectiveLeafTypes @assert precompile(type1,(type2,)) end
end
for type in LinearObjectiveLeafTypes @assert precompile(sparse,(type,)) end
for type in LinearObjectiveLeafTypes @assert precompile(get_numVariables,(type,)) end
for type in LinearObjectiveLeafTypes @assert precompile(get_dependency,(type,)) end
for type in LinearObjectiveLeafTypes @assert precompile(get_gradientSparsity,(type,)) end
for type in LinearObjectiveLeafTypes @assert precompile(get_hessianSparsity,(type,)) end
for type in LinearObjectiveLeafTypes @assert precompile(get_gradientType,(type,)) end
for type in LinearObjectiveLeafTypes @assert precompile(get_hessianType,(type,)) end
for type in LinearObjectiveLeafTypes @assert precompile(firstOrderTaylor,(type,Vector{Float})) end
for type in LinearObjectiveLeafTypes @assert precompile(evaluate,(type,Vector{Float})) end
for type in LinearObjectiveLeafTypes @assert precompile(evaluate_gradient,(type,Vector{Float})) end
for type in LinearObjectiveLeafTypes @assert precompile(evaluate_hessian,(type,Vector{Float})) end
