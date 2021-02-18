# @Author: Massimo De Mauri <massimo>
# @Date:   2019-05-23T18:25:57+02:00
# @Email:  massimo.demauri@gmail.com
# @Filename: LinearConstraintSet.jl
# @Last modified by:   massimo
# @Last modified time: 2021-01-21T11:36:09+01:00
# @License: LGPL-3.0
# @Copyright: {{copyright}}

# named constructor
function LinearConstraintSet(;A::T,loBs::Vector{Float},upBs::Vector{Float})::LinearConstraintSet{T} where T<:Union{Matrix{Float},SpMatrix{Float}}
    return LinearConstraintSet{T}(A,loBs,upBs)
end

# copy functions (Fundamental. Those are used in Branch and Bound)
function Base.copy(cnsSet::LinearConstraintSet)::LinearConstraintSet
    return LinearConstraintSet(cnsSet.A,cnsSet.loBs,cnsSet.upBs)
end

function Base.deepcopy(cnsSet::LinearConstraintSet)::LinearConstraintSet
    return LinearConstraintSet(deepcopy(cnsSet.A),copy(cnsSet.loBs),copy(cnsSet.upBs))
end

# type conversion
function LinearConstraintSet{T}(cnsSet::LinearConstraintSet{T})::LinearConstraintSet{T} where T<:Union{Matrix{Float},SpMatrix{Float}}
    return deepcopy(cnsSet)
end

function LinearConstraintSet{T}(cnsSet::LinearConstraintSet)::LinearConstraintSet{T} where T<:Union{Matrix{Float},SpMatrix{Float}}
    return LinearConstraintSet(T(cnsSet.A),copy(cnsSet.loBs),copy(cnsSet.upBs))
end

function SparseArrays.sparse(cnsSet::LinearConstraintSet)::LinearConstraintSet{SpMatrix{Float}}
    return LinearConstraintSet{SpMatrix{Float}}(cnsSet)
end

# inspect functions (Fundamental. Those are used in Branch and Bound)
function get_numVariables(cnsSet::LinearConstraintSet)::Int
    return size(cnsSet.A,2)
end

function get_size(cnsSet::LinearConstraintSet)::Int
    return size(cnsSet.A,1)
end

function get_bounds(cnsSet::LinearConstraintSet)::Tuple{Vector{Float},Vector{Float}}
    return (cnsSet.loBs,cnsSet.upBs)
end

function Base.getindex(cnsSet::LinearConstraintSet,list::Vector{T})::LinearConstraintSet where T<:Integer
    return LinearConstraintSet(A=cnsSet.A[list,:],
                               loBs=cnsSet.loBs[list],
                               upBs=cnsSet.upBs[list])
end

function Base.getindex(cnsSet::LinearConstraintSet,range::UnitRange)::LinearConstraintSet
    return LinearConstraintSet(A=cnsSet.A[range,:],
                               loBs=cnsSet.loBs[range],
                               upBs=cnsSet.upBs[range])
end

function Base.getindex(cnsSet::LinearConstraintSet,index::T)::LinearConstraintSet where T <: Integer
    return getindex(cnsSet,[index])
end

function get_dependency(cnsSet::LinearConstraintSet)::Vector{Vector{Int}}
    return [get_dependency(cnsSet,index) for index in 1:get_size(cnsSet)]
end

function get_dependency(cnsSet::LinearConstraintSet{Matrix{Float}},index::Int)::Vector{Int}
    return findall(!iszero,cnsSet.A[index,:])
end

function get_dependency(cnsSet::LinearConstraintSet{SpMatrix{Float}},index::Int)::Vector{Int}
    return findnz(cnsSet.A[index,:])[1]
end

function get_firstNZs(cnsSet::LinearConstraintSet,dimension::Int)::Vector{Int}
    @assert 0 < dimension <= 2
    return findFirstNZs(cnsSet.A,dimension)
end

function get_firstNZs(cnsSet::LinearConstraintSet,indices::Vector{Int},dimension::Int)::Vector{Int}
    @assert 0 < dimension <= 2
    if dimension == 1
        return findFirstNZs(cnsSet.A[:,indices],dimension)
    else
        return findFirstNZs(cnsSet.A[indices,:],dimension)
    end
end

function get_lastNZs(cnsSet::LinearConstraintSet,dimension::Int)::Vector{Int}
    @assert 0 < dimension <= 2
    return findLastNZs(cnsSet.A,dimension)
end

function get_lastNZs(cnsSet::LinearConstraintSet,indices::Vector{Int},dimension::Int)::Vector{Int}
    @assert 0 < dimension <= 2
    if dimension == 1
        return findLastNZs(cnsSet.A[:,indices],dimension)
    else
        return findLastNZs(cnsSet.A[indices,:],dimension)
    end
end

function get_jacobianSparsity(cnsSet::LinearConstraintSet)::Tuple{Vector{Int},Vector{Int}}
    if cnsSet isa LinearConstraintSet{Matrix{Float}}
        return findnz(SpMatrix{Float}(cnsSet.A))[1:2]
    else
        return findnz(cnsSet.A)[1:2]
    end
end

function get_hessianSparsity(cnsSet::LinearConstraintSet)::Tuple{Vector{Int},Vector{Int}}
    return (Int[],Int[])
end

function get_jacobianType(cnsSet::LinearConstraintSet{T})::Type where T<:Union{Matrix{Float},SpMatrix{Float}}
    return T
end

function get_hessianType(cnsSet::LinearConstraintSet{T})::Type where T<:Union{Matrix{Float},SpMatrix{Float}}
    return T
end

function linearRelaxation(cnsSet::LinearConstraintSet,linPoint::Vector{Float})::LinearConstraintSet
    return cnsSet
end

# give the value of the constraints in the given point
function evaluate(cnsSet::LinearConstraintSet,point::VirtualVector{Float})::Vector{Float}
    @assert length(point) == get_numVariables(cnsSet)
    return cnsSet.A*point
end

# evaluate the jacobian of the constraints in the given point
function evaluate_jacobian(cnsSet::LinearConstraintSet{T},point::VirtualVector{Float})::T  where T<:Union{Matrix{Float},SpMatrix{Float}}
    @assert length(point) == get_numVariables(cnsSet)
    return copy(cnsSet.A)
end

function evaluate_hessian(cnsSet::LinearConstraintSet{T},point::VirtualVector{Float})::Vector{T} where T<:Union{Matrix{Float},SpMatrix{Float}}
    @assert length(point) == get_numVariables(cnsSet)
    numCnss = get_size(cnsSet)
    numVars = get_numVariables(costraintSet)
    if T == Matrix{Float}
        return [zeros(numVars,numVars) for k in 1:numCnss]
    else
        return [spzeros(numVars,numVars) for k in 1:numCnss]
    end
end


# update functions (Not fundamental. Those are used only in updating the problem)
function update_bounds!(cnsSet::LinearConstraintSet;loBs::Vector{Float}=Float[],upBs::Vector{Float}=Float[])::Nothing
    if length(loBs) > 0
        @assert length(loBs) == length(cnsSet.loBs) == length(cnsSet.upBs)
        @. cnsSet.loBs = loBs
    end
    if length(upBs) > 0
        @assert length(upBs) == length(cnsSet.loBs) == length(cnsSet.upBs)
        @. cnsSet.upBs = upBs
    end
    return
end

function update_bounds!(cnsSet::LinearConstraintSet,indices::Vector{Int};loBs::Vector{Float}=Float[],upBs::Vector{Float}=Float[])::Nothing
    if length(loBs) > 0
        @assert length(loBs) == length(indices)
        @. cnsSet.loBs[indices] = loBs
    end
    if length(upBs) > 0
        @assert length(upBs) == length(indices)
        @. cnsSet.upBs[indices] = upBs
    end
    return
end

function remove_variables!(cnsSet::LinearConstraintSet{T},indices::Union{Vector{Int},UnitRange{Int}})::Nothing where T<:Union{Matrix{Float},SpMatrix{Float}}
    if nnz(cnsSet.A) > 0
        toKeep = filter(x->!(x in indices), collect(1:get_numVariables(cnsSet)))
        cnsSet.A = cnsSet.A[:,toKeep]
    else
        cnsSet.A = cnsSet.A[:,1:end-length(indices)]
    end

    return
end

function fix_variables!(cnsSet::LinearConstraintSet,indices::Union{Vector{Int},UnitRange{Int}},values::Vector{Float};removeFixedVariables::Bool=false)::Nothing
    @assert length(indices) == length(values)

    if nnz(cnsSet.A) > 0
        deltaBounds = -cnsSet.A[:,indices]*values
        cnsSet.loBs += deltaBounds
        cnsSet.upBs += deltaBounds
        if removeFixedVariables
            toKeep = filter(x->!(x in indices), collect(1:get_numVariables(cnsSet)))
            cnsSet.A = cnsSet.A[:,toKeep]
        else
            if cnsSet.A isa Matrix
                cnsSet.A[:,indices] = zeros(size(cnsSet.A,1),length(indices))
            else
                cnsSet.A[:,indices] = spzeros(size(cnsSet.A,1),length(indices))
            end
        end
    else
        if removeFixedVariables
            cnsSet.A = cnsSet.A[:,end-length(indices)+1:end]
        end
    end

    return
end

function insert_variables!(cnsSet::LinearConstraintSet{T},numVariables::Int,insertionPoint::Int)::Nothing where T<:Union{Matrix{Float},SpMatrix{Float}}
    @assert numVariables>=0
    @assert 0<insertionPoint<=get_numVariables(cnsSet)+1
    cnsSet.A = hcat(cnsSet.A[:,1:insertionPoint-1],zeros(size(cnsSet.A,1),numVariables),cnsSet.A[:,insertionPoint:end])
    return
end

function append_variables!(cnsSet::LinearConstraintSet{T},numVariables::Int)::Nothing where T<:Union{Matrix{Float},SpMatrix{Float}}
    insert_variables!(cnsSet,numVariables,get_numVariables(cnsSet)+1)
    return
end

function remove_constraints!(cnsSet::LinearConstraintSet{T},indices::Union{Vector{Int},UnitRange{Int}})::Nothing  where T<:Union{Matrix{Float},SpMatrix{Float}}
    cnsSet.A = cnsSet.A[[i for i in 1:get_size(cnsSet) if !(i in indices)],:]
    deleteat!(cnsSet.loBs,indices)
    deleteat!(cnsSet.upBs,indices)
    return
end


function Base.insert!(cnsSet1::LinearConstraintSet{T1},cnsSet2::LinearConstraintSet{T2},insertionPoint::Int)::Nothing where T1<:Union{Matrix{Float},SpMatrix{Float}} where T2<:Union{Matrix{Float},SpMatrix{Float}}
    cnsSet1.A = vcat(cnsSet1.A[1:insertionPoint-1,:],cnsSet2.A,cnsSet1.A[insertionPoint:end,:])
    splice!(cnsSet1.loBs,insertionPoint:insertionPoint-1,copy(cnsSet2.loBs))
    splice!(cnsSet1.upBs,insertionPoint:insertionPoint-1,copy(cnsSet2.upBs))
    return
end


function Base.permute!(cnsSet::LinearConstraintSet{T},permutation::Vector{Int})::Nothing where T<:Union{Matrix{Float},SpMatrix{Float}}
    cnsSet.A = cnsSet.A[permutation,:]
    cnsSet.loBs = cnsSet.loBs[permutation]
    cnsSet.upBs = cnsSet.upBs[permutation]
    return
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
        A = sparse(Vector{Int}(serial[offset+1:offset+numNZs]),
                   Vector{Int}(serial[offset+numNZs+1:offset+2*numNZs]),
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


function serialSize(cnsSet::LinearConstraintSet{T})::Int where T<:Union{Matrix{Float},SpMatrix{Float}}
    # estimate necessary memory
    numCnss = size(cnsSet.A,1)
    numVars = size(cnsSet.A,2)
    if T == Matrix{Float}
        return 3 + numCnss*numVars + 2*numCnss
    else
        return 3 + 3*nnz(cnsSet.A) + 2*numCnss
    end
end

function serialize(cnsSet::LinearConstraintSet{T})::SerialData where T<:Union{Matrix{Float},SpMatrix{Float}}

    # allocate memory
    serial = SerialData(Vector{Float}(undef,serialSize(cnsSet)))

    # write serialized data
    serialize_in!(serial,cnsSet,offset=0)
    return serial
end


function serialize_in!(serial::SerialData,cnsSet::LinearConstraintSet{T};offset::Int=0)::Int where T<:Union{Matrix{Float},SpMatrix{Float}}

    # check input
    numCnss = size(cnsSet.A,1)
    numVars = size(cnsSet.A,2)

    if T == Matrix{Float}
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




############################# precompilation #########################
for type in LinearConstraintSetLeafTypes @assert precompile(copy,(type,)) end
for type in LinearConstraintSetLeafTypes @assert precompile(deepcopy,(type,)) end
for type1 in LinearConstraintSetLeafTypes
    for type2 in LinearConstraintSetLeafTypes @assert precompile(type1,(type2,)) end
end
for type in LinearConstraintSetLeafTypes @assert precompile(sparse,(type,)) end
for type in LinearConstraintSetLeafTypes @assert precompile(get_numVariables,(type,)) end
for type in LinearConstraintSetLeafTypes @assert precompile(get_size,(type,)) end
for type in LinearConstraintSetLeafTypes @assert precompile(get_bounds,(type,)) end
for type in LinearConstraintSetLeafTypes @assert precompile(get_dependency,(type,)) end
for type in LinearConstraintSetLeafTypes @assert precompile(get_dependency,(type,Int)) end
for type in LinearConstraintSetLeafTypes @assert precompile(get_firstNZs,(type,Int)) end
for type in LinearConstraintSetLeafTypes @assert precompile(get_firstNZs,(type,Vector{Int},Int)) end
for type in LinearConstraintSetLeafTypes @assert precompile(get_lastNZs,(type,Int)) end
for type in LinearConstraintSetLeafTypes @assert precompile(get_lastNZs,(type,Vector{Int},Int)) end
for type in LinearConstraintSetLeafTypes @assert precompile(get_jacobianSparsity,(type,)) end
for type in LinearConstraintSetLeafTypes @assert precompile(get_hessianType,(type,)) end
for type in LinearConstraintSetLeafTypes @assert precompile(get_jacobianType,(type,)) end
for type in LinearConstraintSetLeafTypes @assert precompile(get_hessianSparsity,(type,)) end
for type in LinearConstraintSetLeafTypes @assert precompile(linearRelaxation,(type,Vector{Float})) end
for type in LinearConstraintSetLeafTypes @assert precompile(evaluate,(type,Vector{Float})) end
for type in LinearConstraintSetLeafTypes @assert precompile(evaluate_jacobian,(type,Vector{Float})) end
for type in LinearConstraintSetLeafTypes @assert precompile(evaluate_hessian,(type,Vector{Float})) end
