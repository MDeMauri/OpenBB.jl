# @Author: Massimo De Mauri <massimo>
# @Date:   2019-05-23T18:34:36+02:00
# @Email:  massimo.demauri@gmail.com
# @Filename: QuadraticObjective.jl
# @Last modified by:   massimo
# @Last modified time: 2021-02-12T18:25:38+01:00
# @License: LGPL-3.0
# @Copyright: {{copyright}}

function is_linear(objFun::QuadraticObjective)::Bool
    return all(objFun.Q .== 0.0)
end
function is_quadratic(objFun::QuadraticObjective)::Bool
    return true
end
function is_convex(objFun::QuadraticObjective)::Bool
    if objFun isa QuadraticObjective{Matrix{Float}}
        return  all(eigen(objFun.Q).values .>= -1e-10)
    else
        return  all(eigen(Matrix(objFun.Q)).values .>= -1e-10)
    end
end

# constructors and copy functions (Fundamental. These are used in Branch and Bound)
# named constructors
function QuadraticObjective(;Q::TQ,L::TL,c::Float=0.0,pInvQ::Matrix{Float}=Matrix{Float}(undef,0,0))::QuadraticObjective{TQ,TL} where {TQ<:Union{Matrix{Float},SpMatrix{Float}},TL<:Union{Vector{Float},SpVector{Float}}}
    if isempty(pInvQ)
        if !iszero(Q)
            pInvQ = pinv(Matrix{Float}(Q))
        else
            pInvQ = nothing
        end
    end
    return QuadraticObjective(Q,L,c,pInvQ)
end

# copy functions (Fundamental)
function Base.copy(objFun::QuadraticObjective)::QuadraticObjective
    return QuadraticObjective(objFun.Q,objFun.L,objFun.c,objFun.pInvQ)
end
function Base.deepcopy(objFun::QuadraticObjective)::QuadraticObjective
    if !isnothing(objFun.pInvQ)
        return QuadraticObjective(copy(objFun.Q),copy(objFun.L),objFun.c,copy(objFun.pInvQ))
    else
        return QuadraticObjective(copy(objFun.Q),copy(objFun.L),objFun.c,nothing)
    end
end


# type conversions (Fundamental)
function QuadraticObjective{TQ,TL}(objFun::QuadraticObjective{TQ,TL})::QuadraticObjective{TQ,TL} where {TQ<:Union{Matrix{Float},SpMatrix{Float}},TL<:Union{Vector{Float},SpVector{Float}}}
    return deepcopy(objFun)
end
function QuadraticObjective{TQ,TL}(objFun::QuadraticObjective)::QuadraticObjective{TQ,TL} where {TQ<:Union{Matrix{Float},SpMatrix{Float}},TL<:Union{Vector{Float},SpVector{Float}}}
    if !isnothing(objFun.pInvQ)
        return QuadraticObjective(TQ(objFun.Q),TL(objFun.L),objFun.c,copy(objFun.pInvQ))
    else
        return QuadraticObjective(TQ(objFun.Q),TL(objFun.L),objFun.c,nothing)
    end
end
function QuadraticObjective{TQ,TL}(objFun::LinearObjective)::QuadraticObjective{TQ,TL} where {TQ<:Union{Matrix{Float},SpMatrix{Float}},TL<:Union{Vector{Float},SpVector{Float}}}
    numVars = get_numVariables(objFun)
    if TQ == Matrix{Float}
        return QuadraticObjective(zeros(numVars,numVars),TL(objFun.L),objFun.c,nothing)
    else
        return QuadraticObjective(spzeros(numVars,numVars),TL(objFun.L),objFun.c,nothing)
    end
end


function SparseArrays.sparse(objFun::QuadraticObjective)::QuadraticObjective{SpMatrix{Float},SpVector{Float}}
    return QuadraticObjective{SpMatrix{Float},SpVector{Float}}(objFun)
end


# inspect functions  (Fundamental. These are used in Branch and Bound)
function get_numVariables(objFun::QuadraticObjective)::Int
    return size(objFun.Q,1)
end

function get_dependency(objFun::QuadraticObjective{TQ,TL})::Vector{Int} where {TQ<:Union{Matrix{Float},SpMatrix{Float}},TL<:Union{Vector{Float},SpVector{Float}}}

    if TQ == Matrix{Float}
        Q = sparse(objFun.Q)
    else
        Q = objFun.Q
    end
    if TL == Vector{Float}
        L = sparse(objFun.L)
    else
        L = objFun.L
    end
    dependency = vcat(findnz(Q)[2],findnz(L)[1])
    unique!(sort!(dependency))
    return dependency

end

function get_gradientSparsity(objFun::QuadraticObjective)::Vector{Int}
    numVars = get_numVariables(objFun)
    return get_dependency(objFun)
end

function get_hessianSparsity(objFun::QuadraticObjective)::Tuple{Vector{Int},Vector{Int}}
    return findnz(SpMatrix{Float}(objFun.Q))[1:2]
end

function get_gradientType(constraintSet::QuadraticObjective{<:AbstractMatrix,TL})::Type where TL<:Union{Vector{Float},SpVector{Float}}
    return TL
end

function get_hessianType(constraintSet::QuadraticObjective{TQ})::Type where TQ<:Union{Matrix{Float},SpMatrix{Float}}
    return TQ
end


# evaluate the objFun function in the given point
function evaluate(objFun::QuadraticObjective,point::VirtualVector{Float})::Float
    @assert length(point) == get_numVariables(objFun)
    return .5*point'*objFun.Q*point + objFun.L'*point + objFun.c
end

# evaluate the gradient of the objFun function in the given point
function evaluate_gradient(objFun::QuadraticObjective{TQ,TL},point::VirtualVector{Float})::TL where {TQ<:Union{Matrix{Float},SpMatrix{Float}},TL<:Union{Vector{Float},SpVector{Float}}}
    @assert length(point) == get_numVariables(objFun)
    return objFun.Q*point + objFun.L
end

# evaluate the hessian of the objFun function in the given point
function evaluate_hessian(objFun::QuadraticObjective{TQ,TL},point::VirtualVector{Float})::TQ where {TQ<:Union{Matrix{Float},SpMatrix{Float}},TL<:Union{Vector{Float},SpVector{Float}}}
    @assert length(point) == get_numVariables(objFun)
    return copy(objFun.Q)
end



# update functions (Not fundamental These are used only for problem update)
function insert_variables!(objFun::QuadraticObjective,numVariables::Int,insertionPoint::Int)::Nothing
    @assert numVariables >= 0
    @assert 1<=insertionPoint<=get_numVariables(objFun)+1

    if objFun.Q isa Matrix{Float}
        objFun.Q = vcat(hcat(objFun.Q[1:insertionPoint-1,1:insertionPoint-1],zeros(insertionPoint-1,numVariables),objFun.Q[1:insertionPoint-1,insertionPoint:end]),
                           zeros(numVariables,numVariables+get_numVariables(objFun)),
                           hcat(objFun.Q[insertionPoint:end,1:insertionPoint-1],zeros(get_numVariables(objFun)-insertionPoint+1,numVariables),objFun.Q[insertionPoint:end,insertionPoint:end]))
    else
        objFun.Q = vcat(hcat(objFun.Q[1:insertionPoint-1,1:insertionPoint-1],spzeros(insertionPoint-1,numVariables),objFun.Q[1:insertionPoint-1,insertionPoint:end]),
                           spzeros(numVariables,numVariables+get_numVariables(objFun)),
                           hcat(objFun.Q[insertionPoint:end,1:insertionPoint-1],spzeros(get_numVariables(objFun)-insertionPoint+1,numVariables),objFun.Q[insertionPoint:end,insertionPoint:end]))
    end

    if objFun.L isa Vector{Float}
        objFun.L = vcat(objFun.L[1:insertionPoint-1],zeros(numVariables),objFun.L[insertionPoint:end])
    else
        objFun.L = vcat(objFun.L[1:insertionPoint-1],spzeros(numVariables),objFun.L[insertionPoint:end])
    end

    if !iszero(objFun.Q)
        objFun.pInvQ = pinv(Matrix{Float}(objFun.Q))
    else
        objFun.pInvQ = nothing
    end

    return
end

function append_variables!(objFun::QuadraticObjective,numVariables::Int)::Nothing
    insert_variables!(objFun,numVariables,get_numVariables(objFun)+1)
    return
end


function remove_variables!(objFun::QuadraticObjective,indices::Union{Vector{Int},UnitRange{Int}})::Nothing
    toKeep = filter(x->!(x in indices), collect(1:get_numVariables(objFun)))
    objFun.Q = objFun.Q[toKeep,toKeep]
    objFun.L = objFun.L[toKeep]
    if !iszero(objFun.Q)
        objFun.pInvQ = pinv(Matrix{Float}(objFun.Q))
    else
        objFun.pInvQ = nothing
    end
    return
end

function fix_variables!(objFun::QuadraticObjective,indices::Union{Vector{Int},UnitRange{Int}},values::Vector{Float};removeFixedVariables::Bool=false)::Nothing

    objFun.c += .5*values'*objFun.Q[indices,indices]*values + objFun.L'*values
    toKeep = filter(x->!(x in indices), collect(1:get_numVariables(objFun)))
    if removeFixedVariables
        objFun.L = objFun.L[toKeep] + dropdims(objFun.Q[toKeep,indices]*values,dims=2)
        objFun.Q = objFun.Q[toKeep,toKeep]
    else
        objFun.L[toKeep] += dropdims(objFun.Q[toKeep,indices]*values,dims=2)
        if objFun.L isa Vector
            objFun.L[indices] .= 0.0
        else
            objFun.L[indices] = spzeros(length(indices))
        end
        if objFun.Q isa Matrix
            objFun.Q[:,indices] = zeros(size(objFun.Q,1),length(indices))
            objFun.Q[indices,:] = zeros(length(indices),size(objFun.Q,1))
        else
            objFun.Q[:,indices] = spzeros(size(objFun.Q,1),length(indices))
            objFun.Q[indices,:] = spzeros(length(indices),size(objFun.Q,2))
        end
    end
    if !iszero(objFun.Q)
        objFun.pInvQ = pinv(Matrix{Float}(objFun.Q))
    else
        objFun.pInvQ = nothing
    end
    return
end


import Base.+
function +(objFun1::QuadraticObjective{TQ,TL},objFun2::QuadraticObjective{TQ,TL})::QuadraticObjective{TQ,TL} where {TQ<:Union{Matrix{Float},SpMatrix{Float}},TL<:Union{Vector{Float},SpVector{Float}}}
    @assert get_numVariables(objFun1) == get_numVariables(objFun2)
    return QuadraticObjective(Q=objFun1.Q+objFun2.Q,L=objFun1.L+objFun2.L,c=objFun1.c+objFun2.c)
end

function add!(objFun1::QuadraticObjective{TQ,TL},objFun2::QuadraticObjective{TQ,TL})::Nothing where {TQ<:Union{Matrix{Float},SpMatrix{Float}},TL<:Union{Vector{Float},SpVector{Float}}}
    @assert get_numVariables(objFun1) == get_numVariables(objFun2)
    objFun1.Q += objFun2.Q
    objFun1.L += objFun2.L
    objFun1.c += objFun2.c
    if !iszero(objFun1.Q)
        objFun1.pInvQ = pinv(Matrix{Float}(objFun1.Q))
    else
        objFun1.pInvQ = nothing
    end
    return
end



## Serialization (Not fundamental. Used to store and to send)

function serialSize(objFun::QuadraticObjective{TQ,TL})::Int where {TQ<:Union{Matrix{Float},SpMatrix{Float}},TL<:Union{Vector{Float},SpVector{Float}}}
    # estimate necessary memory
    numVars = size(objFun.Q,1)
    sSize = 4
    if TQ == Matrix{Float}
        sSize += numVars^2
    else
        sSize += 3*nnz(objFun.Q)
    end

    if TL == Vector{Float}
        sSize += numVars
    else
        sSize += 2*nnz(objFun.L)
    end

    return sSize
end


function QuadraticObjective(serial::SerialData;offset::Int=0)::Tuple{QuadraticObjective,Int}
    # header
    numVars = Int(serial[offset+1])
    numNZsQ = Int(serial[offset+2])
    numNZsL = Int(serial[offset+3])
    offset += 3

    if numNZsQ == -1 # dense quadratic term
        # check input
        @assert length(serial) >= offset + numVars^2 + 1
        # reconstruct constraint matrix
        Q = zeros(numVars,numVars)
        for k in 1:numVars
            Q[k,:] = serial[offset+1:offset+numVars]
            offset += numVars
        end
    else # sparse quadratic term
        # check input
        @assert length(serial) >= offset + 3*numNZsQ + 1
        # reconstruct constraint matrix
        Q = sparse(Vector{Int}(serial[offset+1:offset+numNZsQ]),
                   Vector{Int}(serial[offset+numNZsQ+1:offset+2*numNZsQ]),
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
        L = sparsevec(Vector{Int}(serial[offset+1:offset+numNZsL]),serial[offset+numNZsL+1:offset+2*numNZsL],numVars)
        offset += 2*numNZsL
    end

    # constant term
    c = serial[offset+1]
    offset += 1

    # reconstruct the objFun function
    return (QuadraticObjective(Q=Q,L=L,c=c),offset)

end

function serialize(objFun::QuadraticObjective)::SerialData

    # allocate memory
    serial = SerialData(Vector{Float}(undef,serialSize(objFun)))

    # write serialized data
    serialize_in!(serial,objFun,offset=0)
    return serial
end

function serialize_in!(serial::SerialData,objFun::QuadraticObjective{TQ,TL};offset::Int=0)::Int where {TQ<:Union{Matrix{Float},SpMatrix{Float}},TL<:Union{Vector{Float},SpVector{Float}}}
    # check input
    @assert length(serial) >= offset + serialSize(objFun)

    # collect info
    numVars = size(objFun.Q,1)

    # header
    serial[offset+1] = numVars
    if TQ == Matrix{Float}
        serial[offset+2] = -1.0
    else
        serial[offset+2] = nnz(objFun.Q)
    end

    if TL == Vector{Float}
        serial[offset+3] = -1.0
    else
        serial[offset+3] = nnz(objFun.L)
    end
    offset += 3

    if TQ == Matrix{Float} # dense quadratic term
        for k in 1:numVars
            serial[offset+1:offset+numVars] = objFun.Q[k,:]; offset += numVars
        end
    else # sparse quadratic term
        (rows,cols,vals) = findnz(objFun.Q)
        numNZs_ = length(vals)
        serial[offset+1:offset+numNZs_] = rows; offset += numNZs_
        serial[offset+1:offset+numNZs_] = cols; offset += numNZs_
        serial[offset+1:offset+numNZs_] = vals; offset += numNZs_
    end

    if TL == Vector{Float} #dense linear term
        serial[offset+1:offset+numVars] = objFun.L; offset += numVars
    else # sparse linear term
        (inds,vals) = findnz(objFun.L)
        numNZs_ = length(vals)
        serial[offset+1:offset+numNZs_] = inds; offset += numNZs_
        serial[offset+1:offset+numNZs_] = vals; offset += numNZs_
    end

    # constant term
    serial[offset+1] = objFun.c; offset += 1

    return offset
end


for type in QuadraticObjectiveLeafTypes @assert precompile(is_linear,(type,)) end
for type in QuadraticObjectiveLeafTypes @assert precompile(is_quadratic,(type,)) end
for type in QuadraticObjectiveLeafTypes @assert precompile(is_convex,(type,)) end
for type in QuadraticObjectiveLeafTypes @assert precompile(copy,(type,)) end
for type in QuadraticObjectiveLeafTypes @assert precompile(deepcopy,(type,)) end
for type1 in QuadraticObjectiveLeafTypes
    for type2 in QuadraticObjectiveLeafTypes @assert precompile(type1,(type2,)) end
    for type2 in LinearObjectiveLeafTypes @assert precompile(type1,(type2,)) end
end
for type in QuadraticObjectiveLeafTypes @assert precompile(sparse,(type,)) end
for type in QuadraticObjectiveLeafTypes @assert precompile(get_numVariables,(type,)) end
for type in QuadraticObjectiveLeafTypes @assert precompile(get_dependency,(type,)) end
for type in QuadraticObjectiveLeafTypes @assert precompile(get_gradientSparsity,(type,)) end
for type in QuadraticObjectiveLeafTypes @assert precompile(get_hessianSparsity,(type,)) end
for type in QuadraticObjectiveLeafTypes @assert precompile(get_gradientType,(type,)) end
for type in QuadraticObjectiveLeafTypes @assert precompile(get_hessianType,(type,)) end
for type in QuadraticObjectiveLeafTypes @assert precompile(firstOrderTaylor,(type,Vector{Float})) end
for type in QuadraticObjectiveLeafTypes @assert precompile(evaluate,(type,Vector{Float})) end
for type in QuadraticObjectiveLeafTypes @assert precompile(evaluate_gradient,(type,Vector{Float})) end
for type in QuadraticObjectiveLeafTypes @assert precompile(evaluate_hessian,(type,Vector{Float})) end
