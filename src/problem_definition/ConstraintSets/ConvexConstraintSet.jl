# @Author: Massimo De Mauri <massimo>
# @Date:   2020-10-12T13:24:00+02:00
# @Email:  massimo.demauri@gmail.com
# @Filename: ConvexConstraintSet.jl
# @Last modified by:   massimo
# @Last modified time: 2021-02-16T19:21:20+01:00
# @License: LGPL-3.0
# @Copyright: {{copyright}}



# named constructor
function ConvexConstraintSet(;evalVal::Function,evalJcb::Function,evalHes::Function,
                                 loBs::Vector{Float},upBs::Vector{Float},
                                 jcbSparsity::SpMatrix{Bool},hesSparsity::Vector{SpMatrix{Bool}},
                                 typeJcb::Type,typeHes::Type)::ConvexConstraintSet

    @assert typeJcb<:Union{Matrix{Float},SpMatrix{Float}}
    @assert typeHes<:Union{Matrix{Float},SpMatrix{Float}}
    return ConvexConstraintSet{typeHes,typeJcb}(evalVal,evalJcb,evalHes,jcbSparsity,hesSparsity,loBs,upBs)
end


# copy functions (Fundamental. Those are used in Branch and Bound)
function Base.copy(cnsSet::ConvexConstraintSet{Th,Tj})::ConvexConstraintSet{Th,Tj} where {Th<:Union{Matrix{Float},SpMatrix{Float}},Tj<:Union{Matrix{Float},SpMatrix{Float}}}
    return ConvexConstraintSet{Th,Tj}(cnsSet.evalVal,cnsSet.evalJcb,cnsSet.evalHes,
                                         cnsSet.jcbSparsity,cnsSet.hesSparsity,
                                         cnsSet.loBs,cnsSet.upBs)
end



function Base.deepcopy(cnsSet::ConvexConstraintSet{Th,Tj})::ConvexConstraintSet{Th,Tj} where {Th<:Union{Matrix{Float},SpMatrix{Float}},Tj<:Union{Matrix{Float},SpMatrix{Float}}}
    return ConvexConstraintSet{Th,Tj}(deepcopy(cnsSet.evalVal),deepcopy(cnsSet.evalJcb),deepcopy(cnsSet.evalHes),
                                      copy(cnsSet.jcbSparsity),deepcopy(cnsSet.hesSparsity),copy(cnsSet.loBs),copy(cnsSet.upBs))
end

# type conversion
function ConvexConstraintSet{Th,Tj}(cnsSet::ConvexConstraintSet{Th,Tj})::ConvexConstraintSet{Th,Tj} where {Th<:Union{Matrix{Float},SpMatrix{Float}},Tj<:Union{Matrix{Float},SpMatrix{Float}}}
    return deepcopy(cnsSet)
end

function ConvexConstraintSet{Th,Tj}(cnsSet::ConvexConstraintSet{Th,Tj2})::ConvexConstraintSet{Th,Tj} where {Th<:Union{Matrix{Float},SpMatrix{Float}},Tj<:Union{Matrix{Float},SpMatrix{Float}},Tj2<:Union{Matrix{Float},SpMatrix{Float}}}
    oldEvalVal = cnsSet.evalVal
    oldEvalJcb = cnsSet.evalJcb
    oldEvalHes = cnsSet.evalHes
    return ConvexConstraintSet{Th,Tj}(oldEvalVal,x::VirtualVector{Float}->Tj(oldEvalJcb(x)),cnsSet.evalHes,
                                         deepcopy(cnsSet.jcbSparsity),deepcopy(cnsSet.hesSparsity),
                                         copy(cnsSet.loBs),copy(cnsSet.upBs))
end
function ConvexConstraintSet{Th,Tj}(cnsSet::ConvexConstraintSet{Th2,Tj})::ConvexConstraintSet{Th,Tj} where {Th<:Union{Matrix{Float},SpMatrix{Float}},Tj<:Union{Matrix{Float},SpMatrix{Float}},Th2<:Union{Matrix{Float},SpMatrix{Float}}}
    oldEvalVal = cnsSet.evalVal
    oldEvalJcb = cnsSet.evalJcb
    oldEvalHes = cnsSet.evalHes
    return ConvexConstraintSet{Th,Tj}(oldEvalVal,oldEvalJcb,x::VirtualVector{Float}->Th.(oldEvalHes(x)),
                                         deepcopy(cnsSet.jcbSparsity),deepcopy(cnsSet.hesSparsity),
                                         copy(cnsSet.loBs),copy(cnsSet.upBs))
end
function ConvexConstraintSet{Th,Tj}(cnsSet::ConvexConstraintSet)::ConvexConstraintSet{Th,Tj} where {Th<:Union{Matrix{Float},SpMatrix{Float}},Tj<:Union{Matrix{Float},SpMatrix{Float}}}
    oldEvalVal = cnsSet.evalVal
    oldEvalJcb = cnsSet.evalJcb
    oldEvalHes = cnsSet.evalHes
    return ConvexConstraintSet{Th,Tj}(oldEvalVal,x::VirtualVector{Float}->Tj(oldEvalJcb(x)),x::VirtualVector{Float}->Th.(oldEvalHes(x)),
                                         deepcopy(cnsSet.jcbSparsity),deepcopy(cnsSet.hesSparsity),
                                         copy(cnsSet.loBs),copy(cnsSet.upBs))
end

function ConvexConstraintSet{Th,Tj}(cnsSet::LinearConstraintSet)::ConvexConstraintSet{Th,Tj} where {Th<:Union{Matrix{Float},SpMatrix{Float}},Tj<:Union{Matrix{Float},SpMatrix{Float}}}

    # collect info
    numVars = get_numVariables(cnsSet)
    numCnss = get_size(cnsSet)

    # build the evaluation function
    evalVal = x::VirtualVector{Float}->(A*x)::Vector{Float}

    # build the jacobian and hessian functions
    if Tj == Matrix{Float}
        A = Matrix{Float}(cnsSet.A)
        evalJcb = x::VirtualVector{Float}->(A)::Matrix{Float}
    else
        A = SpMatrix{Float}(cnsSet.A)
        evalJcb = x::VirtualVector{Float}->(A)::SpMatrix{Float}
    end
    if Th == Matrix{Float}
        evalHes = x::VirtualVector{Float}->[zeros(numVars,numVars) for k in 1:numCnss]::Vector{Matrix{Float}}
    else
        evalHes = x::VirtualVector{Float}->[spzeros(numVars,numVars) for k in 1:numCnss]::Vector{SpMatrix{Float}}
    end

    # build jacobian and hessian sparsity
    jcbSparsity=sparsityMatrix(cnsSet.A)
    hesSparsity=[spzeros(Bool,numVars,numVars) for k in 1:numCnss]

    return ConvexConstraintSet{Th,Tj}(evalVal,evalJcb,evalHes,jcbSparsity,hesSparsity,cnsSet.loBs,cnsSet.upBs)
end


function SparseArrays.sparse(cnsSet::ConvexConstraintSet)::ConvexConstraintSet{SpMatrix{Float},SpMatrix{Float}}
    return ConvexConstraintSet{SpMatrix{Float},SpMatrix{Float}}(cnsSet)
end


# inspect functions (Fundamental. Those are used in Branch and Bound)
function get_numVariables(cnsSet::ConvexConstraintSet)::Int
    return size(cnsSet.jcbSparsity,2)
end

function get_size(cnsSet::ConvexConstraintSet)::Int
    return size(cnsSet.jcbSparsity,1)
end

function get_bounds(cnsSet::ConvexConstraintSet)::Tuple{Vector{Float},Vector{Float}}
    return (cnsSet.loBs,cnsSet.upBs)
end

function get_dependency(cnsSet::ConvexConstraintSet)::Vector{Vector{Int}}
    return [get_dependency(cnsSet,index) for index in 1:get_size(cnsSet)]
end


function get_dependency(cnsSet::ConvexConstraintSet,index::Int)::Vector{Int}
    return findnz(cnsSet.jcbSparsity[index,:])[1]
end

function get_firstNZs(cnsSet::ConvexConstraintSet,dimension::Int)::Vector{Int}
    @assert 1 <= dimension <= 2
    return findFirstNZs(cnsSet.jcbSparsity,dimension)
end

function get_firstNZs(cnsSet::ConvexConstraintSet,indices::Vector{Int},dimension::Int)::Vector{Int}
    @assert 1 <= dimension <= 2
    if dimension == 1
        return findFirstNZs(cnsSet.jcbSparsity[:,indices],1)
    else
        return findFirstNZs(cnsSet.jcbSparsity[indices,:],2)
    end
end

function get_lastNZs(cnsSet::ConvexConstraintSet,dimension::Int)::Vector{Int}
    @assert 1 <= dimension <= 2
    return findLastNZs(cnsSet.jcbSparsity,dimension)
end

function get_lastNZs(cnsSet::ConvexConstraintSet,indices::Vector{Int},dimension::Int)::Vector{Int}
    @assert 1 <= dimension <= 2
    if dimension == 1
        return findLastNZs(cnsSet.jcbSparsity[:,indices],dimension)
    else
        return findLastNZs(cnsSet.jcbSparsity[indices,:],dimension)
    end
end

function get_jacobianSparsity(cnsSet::ConvexConstraintSet)::Tuple{Vector{Int},Vector{Int}}
    return findnz(cnsSet.jcbSparsity)[1:2]
end

function get_hessianSparsity(cnsSet::ConvexConstraintSet)::Vector{Tuple{Vector{Int},Vector{Int}}}
    return [findnz(cnsSet.hesSparsity[k])[1:2] for k in 1:get_size(cnsSet)]
end

function get_jacobianType(cnsSet::ConvexConstraintSet{<:AbstractMatrix,Tj})::Type where Tj<:Union{Matrix{Float},SpMatrix{Float}}
    return Tj
end

function get_hessianType(cnsSet::ConvexConstraintSet{Th})::Type where Th<:Union{Matrix{Float},SpMatrix{Float}}
    return Th
end

function linearRelaxation(cnsSet::ConvexConstraintSet,linPoint::Vector{Float})::LinearConstraintSet
    F = cnsSet.evalVal(linPoint)
    J = cnsSet.evalJcb(linPoint)
    deltaBs = -F + J*linPoint
    return LinearConstraintSet(A=J,loBs=cnsSet.loBs .+ deltaBs,upBs=cnsSet.upBs .+ deltaBs)
end

# give the value of the constraints in the given point
function evaluate(cnsSet::ConvexConstraintSet{Th,Tj},point::VirtualVector{Float})::Vector{Float}where {Th<:Union{Matrix{Float},SpMatrix{Float}},Tj<:Union{Matrix{Float},SpMatrix{Float}}}
    @assert length(point) == get_numVariables(cnsSet)
    return cnsSet.evalVal(point)
end

# evaluate the jacobian of the constraints in the given point
function evaluate_jacobian(cnsSet::ConvexConstraintSet{Th,Tj},point::VirtualVector{Float})::Tj where {Th<:Union{Matrix{Float},SpMatrix{Float}},Tj<:Union{Matrix{Float},SpMatrix{Float}}}
    @assert length(point) == get_numVariables(cnsSet)
    return cnsSet.evalJcb(point)
end

# evaluate the hessian of each constraint separately
function evaluate_hessian(cnsSet::ConvexConstraintSet{Th,Tj},point::VirtualVector{Float})::Vector{Th} where {Th<:Union{Matrix{Float},SpMatrix{Float}},Tj<:Union{Matrix{Float},SpMatrix{Float}}}
    @assert length(point) == get_numVariables(cnsSet)
    return cnsSet.evalHes(point)
end

# update functions (Not fundamental. Those are used only in updating the problem)
function update_bounds!(cnsSet::ConvexConstraintSet;loBs::Vector{Float}=Float[],upBs::Vector{Float}=Float[])::Nothing
    if length(loBs) > 0
        @assert length(loBs) == length(cnsSet.upBs)
        @. cnsSet.loBs = loBs
    end
    if length(upBs) > 0
        @assert length(upBs) == length(cnsSet.upBs)
        @. cnsSet.upBs = upBs
    end
    return
end

function update_bounds!(cnsSet::ConvexConstraintSet,indices::Vector{Int};loBs::Vector{Float}=Float[],upBs::Vector{Float}=Float[])::Nothing
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

function fix_variables!(cnsSet::ConvexConstraintSet{Th,Tj},indices::Union{Vector{Int},UnitRange{Int}},values::Vector{Float};removeFixedVariables::Bool=false)::Nothing where {Th<:Union{Matrix{Float},SpMatrix{Float}},Tj<:Union{Matrix{Float},SpMatrix{Float}}}

    # collect info
    oldNumVars = get_numVariables(cnsSet)
    toKeep = filter(x->!(x in indices), collect(1:oldNumVars))
    numCnss = get_size(cnsSet)

    # update evaluation function
    oldEvalVal = cnsSet.evalVal
    x_ = Vector{Float}(undef,oldNumVars); x_[indices] = values
    if removeFixedVariables
        function newEvalVal1(x::VirtualVector{Float})::Vector{Float}
            @. x_[toKeep] = x
            return oldEvalVal(x_)
        end
        cnsSet.evalVal = newEvalVal1
    else
        function newEvalVal2(x::VirtualVector{Float})::Vector{Float}
            @. x_[toKeep] = x[toKeep]
            return oldEvalVal(x_)
        end
        cnsSet.evalVal = newEvalVal2
    end


    # update jacobian
    if Tj == Matrix{Float}
        oldEvalJcb = cnsSet.evalJcb
        x_ = Vector{Float}(undef,oldNumVars); x_[indices] = values
        if removeFixedVariables
            function newEvalJcbD1(x::VirtualVector{Float})::Matrix{Float}
                @. x_[toKeep] = x
                return oldEvalJcb(x_)[:,toKeep]
            end
            cnsSet.evalJcb = newEvalJcbD1
        else
            function newEvalJcbD2(x::VirtualVector{Float})::Matrix{Float}
                @. x_[toKeep] = x[toKeep]
                jacobian = oldEvalJcb(x_)
                jacobian[:,indices] = zeros(size(jacobian,1),length(indices))
                return jacobian
            end
            cnsSet.evalJcb = newEvalJcbD2
        end
    else
        oldEvalJcb = cnsSet.evalJcb
        x_ = Vector{Float}(undef,oldNumVars); x_[indices] = values
        if removeFixedVariables
            function newEvalJcbS1(x::VirtualVector{Float})::SpMatrix{Float}
                @. x_[toKeep] = x
                return oldEvalJcb(x_)[:,toKeep]
            end
            cnsSet.evalJcb = newEvalJcbS1
        else
            function newEvalJcbS2(x::VirtualVector{Float})::SpMatrix{Float}
                @. x_[toKeep] = x[toKeep]
                jacobian = oldEvalJcb(x_)
                jacobian[:,indices] = spzeros(size(jacobian,1),length(indices))
                return jacobian
            end
            cnsSet.evalJcb = newEvalJcbS2
        end
    end


    # update hessian
    if Th == Matrix{Float}
        oldEvalHes = cnsSet.evalHes
        x_ = Vector{Float}(undef,oldNumVars); x_[indices] = values
        if removeFixedVariables
            function newEvalHesD1(x::VirtualVector{Float})::Vector{Matrix{Float}}
                @. x_[toKeep] = x
                hessian_ = oldEvalHes(x_)
                return getindex.(hessian_,[toKeep],[toKeep])
            end
            cnsSet.evalHes = newEvalHesD1
        else
            function newEvalHesD2(x::VirtualVector{Float})::Vector{Matrix{Float}}
                @. x_[toKeep] = x[toKeep]
                hessian_ = oldEvalHes(x_)
                for k in 1:length(hessian_)
                    hessian_[k][:,indices] = zeros(size(hessian_,1),length(indices))
                    hessian_[k][indices,:] = zeros(length(indices),size(hessian_,2))
                end
                return hessian_
            end
            cnsSet.evalHes = newEvalHesD2
        end
    else
        oldEvalHes = cnsSet.evalHes
        x_ = Vector{Float}(undef,oldNumVars); x_[indices] = values
        if removeFixedVariables
            function newEvalHesS1(x::VirtualVector{Float})::Vector{SpMatrix{Float}}
                @. x_[toKeep] = x
                hessian_ = oldEvalHes(x_)
                return getindex.(hessian_,[toKeep],[toKeep])
            end
            cnsSet.evalHes = newEvalHesS1
        else
            function newEvalHesS2(x::VirtualVector{Float})::Vector{SpMatrix{Float}}
                @. x_[toKeep] = x[toKeep]
                hessian_ = oldEvalHes(x_)
                for k in 1:length(hessian_)
                    hessian_[k][:,indices] = spzeros(size(hessian_,1),length(indices))
                    hessian_[k][indices,:] = spzeros(length(indices),size(hessian_,2))
                end
                return hessian_
            end
            cnsSet.evalHes = newEvalHesS2
        end
    end

    # update sparsity matrices
    if removeFixedVariables
        cnsSet.jcbSparsity = cnsSet.jcbSparsity[:,toKeep]
    else
        cnsSet.jcbSparsity[:,indices] = spzeros(Bool,size(cnsSet.jcbSparsity,1),length(indices))
    end
    if removeFixedVariables
        for k in 1:numCnss
            cnsSet.hesSparsity[k] = cnsSet.hesSparsity[k][toKeep,toKeep]
        end
    else
        for k in 1:numCnss
            cnsSet.hesSparsity[k][:,indices] = spzeros(Bool,size(cnsSet.hesSparsity[k],1),length(indices))
            cnsSet.hesSparsity[k][indices,:] = spzeros(Bool,length(indices),size(cnsSet.hesSparsity[k],2))
        end
    end

    return
end


function remove_variables!(cnsSet::ConvexConstraintSet{Th,Tj},indices::Union{Vector{Int},UnitRange{Int}})::Nothing where {Th<:Union{Matrix{Float},SpMatrix{Float}},Tj<:Union{Matrix{Float},SpMatrix{Float}}}

    # check correctness of input
    if !iszero(cnsSet.jcbSparsity[:,indices])
        error("The constraint set depends on the variables to be removed: unsafe to remove them")
        return
    end

    # collect info
    oldNumVars = get_numVariables(cnsSet)
    toKeep = filter(x->!(x in indices), collect(1:oldNumVars))
    numCnss = get_size(cnsSet)

    # update evaluation function
    oldEvalVal = cnsSet.evalVal
    x_ = Vector{Float}(undef,oldNumVars)
    function newEvalVal(x::VirtualVector{Float})::Vector{Float}
        @. x_[toKeep] = x
        return oldEvalVal(x_)
    end
    cnsSet.evalVal = newEvalVal

    # update jacobian
    if Tj == Matrix{Float}
        oldEvalJcb = cnsSet.evalJcb
        x_ = Vector{Float}(undef,oldNumVars)
        function newEvalJcbD(x::VirtualVector{Float})::Matrix{Float}
            @. x_[toKeep] = x
            return oldEvalJcb(x_)[:,toKeep]
        end
        cnsSet.evalJcb = newEvalJcbD
    else
        oldEvalJcb = cnsSet.evalJcb
        x_ = Vector{Float}(undef,oldNumVars)
        function newEvalJcbS(x::VirtualVector{Float})::SpMatrix{Float}
            @. x_[toKeep] = x
            return oldEvalJcb(x_)[:,toKeep]
        end
        cnsSet.evalJcb = newEvalJcbS
    end


    # update hessian
    if Th == Matrix{Float}
        oldEvalHes = cnsSet.evalHes
        x_ = Vector{Float}(undef,oldNumVars)
        function newEvalHesD(x::VirtualVector{Float})::Vector{Matrix{Float}}
            @. x_[toKeep] = x
            hessian_ = oldEvalHes(x_)
            return getindex.(hessian_,[toKeep],[toKeep])
        end
        cnsSet.evalHes = newEvalHesD
    else
        oldEvalHes = cnsSet.evalHes
        x_ = Vector{Float}(undef,oldNumVars)
        function newEvalHesS(x::VirtualVector{Float})::Vector{SpMatrix{Float}}
            @. x_[toKeep] = x
            hessian_ = oldEvalHes(x_)
            return getindex.(hessian_,[toKeep],[toKeep])
        end
        cnsSet.evalHes = newEvalHesS
    end

    # update sparsity matrices
    cnsSet.jcbSparsity = cnsSet.jcbSparsity[:,toKeep]
    for k in 1:numCnss
        cnsSet.hesSparsity[k] = cnsSet.hesSparsity[k][toKeep,toKeep]
    end
end



function insert_variables!(cnsSet::ConvexConstraintSet{Th,Tj},numNewVariables::Int,insertionPoint::Int)::Nothing where {Th<:Union{Matrix{Float},SpMatrix{Float}},Tj<:Union{Matrix{Float},SpMatrix{Float}}}
    @assert numNewVariables>=0
    @assert 1<=insertionPoint<=get_numVariables(cnsSet)+1

    # get info
    numCnss = get_size(cnsSet)
    oldNumVars = get_numVariables(cnsSet)
    newNumVars = oldNumVars + numNewVariables
    indices = vcat(collect(1:insertionPoint-1),collect(insertionPoint+numNewVariables:get_numVariables(cnsSet)+numNewVariables))

    # update evaluation function
    oldEvalVal = cnsSet.evalVal
    function newEvalVal(x::VirtualVector{Float})::Vector{Float}
        return oldEvalVal(view(x,indices))
    end
    cnsSet.evalVal = newEvalVal

    # update jacobian
    if Tj == Matrix{Float}
        oldEvalJcb = cnsSet.evalJcb
        function newEvalJcbD(x::VirtualVector{Float})::Matrix{Float}
            jacobian = oldEvalJcb(view(x,indices))
            return hcat(jacobian[:,1:insertionPoint-1],zeros(numCnss,numNewVariables),jacobian[:,insertionPoint:end])
        end
        cnsSet.evalJcb = newEvalJcbD
    else
        oldEvalJcb = cnsSet.evalJcb
        function newEvalJcbS(x::VirtualVector{Float})::SpMatrix{Float}
            jacobian = oldEvalJcb(view(x,indices))
            return hcat(jacobian[:,1:insertionPoint-1],spzeros(numCnss,numNewVariables),jacobian[:,insertionPoint:end])
        end
        cnsSet.evalJcb = newEvalJcbS
    end

    # update hessian
    if Th == Matrix{Float}
        oldEvalHes = cnsSet.evalHes
        function newEvalHesD(x::VirtualVector{Float})::Vector{Matrix{Float}}
            hessian_ = [zeros(newNumVars,newNumVars) for k in 1:numCnss]
            setindex!.(hessian_,oldEvalHes(view(x,indices)),[indices],[indices])
            return hessian_
        end
        cnsSet.evalHes = newEvalHesD
    else
        oldEvalHes = cnsSet.evalHes
        function newEvalHesS(x::VirtualVector{Float})::Vector{SpMatrix{Float}}
            hessian_ = [spzeros(newNumVars,newNumVars) for k in 1:numCnss]
            setindex!.(hessian_,oldEvalHes(view(x,indices)),[indices],[indices])
            return hessian_
        end
        cnsSet.evalHes = newEvalHesS
    end

    # update sparsity matrices
    newSparsity = spzeros(Bool,numCnss,newNumVars)
    newSparsity[:,indices] = cnsSet.jcbSparsity
    cnsSet.jcbSparsity = newSparsity
    for k in 1:numCnss
        newSparsity = spzeros(Bool,newNumVars,newNumVars)
        newSparsity[indices,indices] = cnsSet.hesSparsity[k]
        cnsSet.hesSparsity[k] = newSparsity
    end

    return
end

function append_variables!(cnsSet::ConvexConstraintSet,numVariables::Int)::Nothing
    insert_variables!(cnsSet,numVariables,get_numVariables(cnsSet)+1)
    return
end

function remove_constraints!(cnsSet::ConvexConstraintSet{Th,Tj},indices::Union{Vector{Int},UnitRange{Int}})::Nothing where {Th<:Union{Matrix{Float},SpMatrix{Float}},Tj<:Union{Matrix{Float},SpMatrix{Float}}}

    # collect info
    toKeep = [i for i in 1:get_size(cnsSet) if !(i in indices)]

    # update value function
    oldEvalVal = cnsSet.evalVal
    function newEvalVal(x::VirtualVector{Float})::Vector{Float}
        return oldEvalVal(x)[toKeep]
    end
    cnsSet.evalVal = newEvalVal

    # update jacobian function
    oldEvalJcb = cnsSet.evalJcb
    function newEvalJcb(x::VirtualVector{Float})::Tj
        return oldEvalJcb(x)[toKeep,:]
    end
    cnsSet.evalJcb = newEvalJcb

    # update hessian function
    oldEvalHes = cnsSet.evalHes
    function newEvalHes(x::VirtualVector{Float})::Vector{Th}
        return oldEvalHes(x)[toKeep]
    end
    cnsSet.evalHes = newEvalHes

    # update the rest
    deleteat!(cnsSet.loBs,indices)
    deleteat!(cnsSet.upBs,indices)
    cnsSet.jcbSparsity = cnsSet.jcbSparsity[toKeep,:]
    deleteat!(cnsSet.hesSparsity,indices)

    return
end

function Base.insert!(cnsSet1::ConvexConstraintSet{Th,Tj},cnsSet2::ConvexConstraintSet{Th,Tj},insertionPoint::Int)::Nothing where {Th<:Union{Matrix{Float},SpMatrix{Float}},Tj<:Union{Matrix{Float},SpMatrix{Float}}}

    # ensure size consistency
    @assert get_numVariables(cnsSet1) == get_numVariables(cnsSet2)

    # update value function
    (oldEvalVal1,oldEvalVal2) = (cnsSet1.evalVal,cnsSet2.evalVal)
    function newEvalVal(x::VirtualVector{Float})::Vector{Float}
        (F1,F2) = (oldEvalVal1(x),oldEvalVal2(x))
        return vcat(F1[1:insertionPoint-1],F2,F1[insertionPoint:end])
    end
    cnsSet1.evalVal = newEvalVal

    # update jacobian
    (oldEvalJcb1,oldEvalJcb2) =(cnsSet1.evalJcb,cnsSet2.evalJcb)
    function newEvalJcb(x::VirtualVector{Float})::Tj
        (J1,J2) = (oldEvalJcb1(x),oldEvalJcb2(x))
        return vcat(J1[1:insertionPoint-1,:],J2,J1[insertionPoint:end,:])
    end
    cnsSet1.evalJcb = newEvalJcb

    # update hessian
    (oldEvalHes1,oldEvalHes2) = (cnsSet1.evalHes,cnsSet2.evalHes)
    function newEvalHes(x::VirtualVector{Float})::Vector{Th}
        (H1,H2) = (oldEvalHes1(x),oldEvalHes2(x))
        return vcat(H1[1:insertionPoint-1],H2,H1[insertionPoint:end])
    end
    cnsSet1.evalHes = newEvalHes

    # update jacobian sparsity
    (SM1,SM2) = (cnsSet1.jcbSparsity,cnsSet2.jcbSparsity)
    cnsSet1.jcbSparsity = vcat(SM1[1:insertionPoint-1,:],SM2,SM1[insertionPoint:end,:])

    # stack hessians sparsities and bounds
    splice!(cnsSet1.hesSparsity,insertionPoint:insertionPoint-1,deepcopy(cnsSet2.hesSparsity))
    splice!(cnsSet1.loBs,insertionPoint:insertionPoint-1,copy(cnsSet2.loBs))
    splice!(cnsSet1.upBs,insertionPoint:insertionPoint-1,copy(cnsSet2.upBs))

    return
end


function Base.permute!(cnsSet::ConvexConstraintSet{Th,Tj},permutation::Vector{Int})::Nothing where {Th<:Union{Matrix{Float},SpMatrix{Float}},Tj<:Union{Matrix{Float},SpMatrix{Float}}}


    # update evaluation func
    oldEvalVal = cnsSet.evalVal
    function newEvalVal(x::VirtualVector{Float})::Vector{Float}
        return oldEvalVal(x)[permutation]
    end
    cnsSet.evalVal = newEvalVal

    # update jacobian
    oldEvalJcb = cnsSet.evalJcb
    function newEvalJcb(x::VirtualVector{Float})::Tj
        return oldEvalJcb(x)[permutation,:]
    end
    cnsSet.evalJcb = newEvalJcb

    # update hessian
    oldEvalHes = cnsSet.evalHes
    function newEvalHes(x::VirtualVector{Float})::Vector{Th}
        return oldEvalHes(x)[permutation]
    end
    cnsSet.evalHes = newEvalHes

    # update sparsities and bounds
    cnsSet.jcbSparsity = cnsSet.jcbSparsity[permutation,:]
    cnsSet.hesSparsity = cnsSet.hesSparsity[permutation]
    cnsSet.loBs = cnsSet.loBs[permutation]
    cnsSet.upBs = cnsSet.upBs[permutation]

    return
end


############################# precompilation #########################
for type in ConvexConstraintSetLeafTypes @assert precompile(copy,(type,)) end
for type in ConvexConstraintSetLeafTypes @assert precompile(deepcopy,(type,)) end
for type1 in ConvexConstraintSetLeafTypes
    for type2 in ConvexConstraintSetLeafTypes @assert precompile(type1,(type2,)) end
    for type2 in LinearConstraintSetLeafTypes @assert precompile(type1,(type2,)) end
end
for type in ConvexConstraintSetLeafTypes @assert precompile(sparse,(type,)) end
for type in ConvexConstraintSetLeafTypes @assert precompile(get_numVariables,(type,)) end
for type in ConvexConstraintSetLeafTypes @assert precompile(get_size,(type,)) end
for type in ConvexConstraintSetLeafTypes @assert precompile(get_bounds,(type,)) end
for type in ConvexConstraintSetLeafTypes @assert precompile(get_dependency,(type,)) end
for type in ConvexConstraintSetLeafTypes @assert precompile(get_dependency,(type,Int)) end
for type in ConvexConstraintSetLeafTypes @assert precompile(get_firstNZs,(type,Int)) end
for type in ConvexConstraintSetLeafTypes @assert precompile(get_firstNZs,(type,Vector{Int},Int)) end
for type in ConvexConstraintSetLeafTypes @assert precompile(get_lastNZs,(type,Int)) end
for type in ConvexConstraintSetLeafTypes @assert precompile(get_lastNZs,(type,Vector{Int},Int)) end
for type in ConvexConstraintSetLeafTypes @assert precompile(get_jacobianSparsity,(type,)) end
for type in ConvexConstraintSetLeafTypes @assert precompile(get_hessianSparsity,(type,)) end
for type in ConvexConstraintSetLeafTypes @assert precompile(get_jacobianType,(type,)) end
for type in ConvexConstraintSetLeafTypes @assert precompile(get_hessianType,(type,)) end
for type in ConvexConstraintSetLeafTypes @assert precompile(linearRelaxation,(type,Vector{Float})) end
for type in ConvexConstraintSetLeafTypes @assert precompile(evaluate,(type,Vector{Float})) end
for type in ConvexConstraintSetLeafTypes @assert precompile(evaluate_jacobian,(type,Vector{Float})) end
for type in ConvexConstraintSetLeafTypes @assert precompile(evaluate_hessian,(type,Vector{Float})) end
