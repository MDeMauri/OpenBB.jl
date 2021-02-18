# @Author: Massimo De Mauri <massimo>
# @Date:   2020-10-23T10:36:10+02:00
# @Email:  massimo.demauri@gmail.com
# @Filename: ConvexObjective.jl
# @Last modified by:   massimo
# @Last modified time: 2021-01-19T14:34:05+01:00
# @License: LGPL-3.0
# @Copyright: {{copyright}}

function is_linear(objective::ConvexObjective)::Bool
    return iszero(objective.hesSparsity)
end
function is_quadratic(objective::ConvexObjective)::Bool
    return false
end
function is_convex(objective::ConvexObjective)::Bool
    return true
end

# named constructor
function ConvexObjective(;evalVal::Function,evalGrd::Function,evalHes::Function,
                          grdSparsity::SpVector{Bool},hesSparsity::SpMatrix{Bool},typeGrd::Type,typeHes::Type
                           )::ConvexObjective

    @assert typeGrd<:Union{Vector{Float},SpVector{Float}}
    @assert typeHes<:Union{Matrix{Float},SpMatrix{Float}}
    return ConvexObjective{typeHes,typeGrd}(evalVal,evalGrd,evalHes,grdSparsity,hesSparsity)
end

# copy functions (Fundamental)
function Base.copy(objective::ConvexObjective{Th,Tg})::ConvexObjective{Th,Tg} where {Th<:Union{Matrix{Float},SpMatrix{Float}},Tg<:Union{Vector{Float},SpVector{Float}}}
    return ConvexObjective(objective.evalVal,objective.evalGrd,objective.evalHes,objective.grdSparsity,objective.hesSparsity)
end
function Base.deepcopy(objective::ConvexObjective{Th,Tg})::ConvexObjective{Th,Tg} where {Th<:Union{Matrix{Float},SpMatrix{Float}},Tg<:Union{Vector{Float},SpVector{Float}}}
    return ConvexObjective{Th,Tg}(deepcopy(objective.evalVal),deepcopy(objective.evalGrd),deepcopy(objective.evalHes),copy(objective.grdSparsity),copy(objective.hesSparsity))
end


# type conversions
function ConvexObjective{Th,Tg}(objective::ConvexObjective{Th,Tg})::ConvexObjective{Th,Tg} where {Th<:Union{Matrix{Float},SpMatrix{Float}},Tg<:Union{Vector{Float},SpVector{Float}}}
    return objective
end
function ConvexObjective{Th,Tg}(objective::ConvexObjective{Th1,Tg})::ConvexObjective{Th,Tg} where {Th<:Union{Matrix{Float},SpMatrix{Float}},Tg<:Union{Vector{Float},SpVector{Float}},Th1<:Union{Matrix{Float},SpMatrix{Float}}}
    oldEvalHes = objective.evalHes
    return ConvexObjective{Th,Tg}(objective.evalVal,objective.evalGrd,x::Vector{Float}->Th(oldEvalHes(x)),copy(objective.grdSparsity),copy(objective.hesSparsity))
end
function ConvexObjective{Th,Tg}(objective::ConvexObjective{Th,Tg1})::ConvexObjective{Th,Tg} where {Th<:Union{Matrix{Float},SpMatrix{Float}},Tg<:Union{Vector{Float},SpVector{Float}},Tg1<:Union{Vector{Float},SpVector{Float}}}
    oldEvalGrd = objective.evalGrd
    return ConvexObjective{Th,Tg}(objective.evalVal,x::Vector{Float}->Tg(oldEvalGrd(x)),objective.evalHes,copy(objective.grdSparsity),copy(objective.hesSparsity))
end
function ConvexObjective{Th,Tg}(objective::ConvexObjective)::ConvexObjective{Th,Tg} where {Th<:Union{Matrix{Float},SpMatrix{Float}},Tg<:Union{Vector{Float},SpVector{Float}}}
    oldEvalHes = objective.evalHes
    oldEvalGrd = objective.evalGrd
    return ConvexObjective{Th,Tg}(objective.evalVal,x::Vector{Float}->Tg(oldEvalGrd(x)),x::Vector{Float}->Th(oldEvalHes(x)),copy(objective.grdSparsity),copy(objective.hesSparsity))
end

function ConvexObjective{Th,Tg}(objective::QuadraticObjective)::ConvexObjective{Th,Tg} where {Th<:Union{Matrix{Float},SpMatrix{Float}},Tg<:Union{Vector{Float},SpVector{Float}}}

    # extract the data from the structure to obtain a correct closure
    Q = Th(objective.Q)
    L = Tg(objective.L)
    c = objective.c

    # define value function, gradient and hessian
    newEvalVal = x::Vector{Float}->(0.5*x'*Q*x + L'*x + c)::Float
    if Tg == SpVector{Float}
        newEvalGrd = x::Vector{Float}->(TQ*x + L)::Tg
    else
        newEvalGrd = x::Vector{Float}->(Q*x + L)::Tg
    end
    newEvalHes = x::Vector{Float}->Q::Th

    return ConvexObjective{Th,Tg}(newEvalVal,newEvalGrd,newEvalHes,sparsityVector(objective.L),sparsityMatrix(objective.Q))
end


function ConvexObjective{Th,Tg}(objective::LinearObjective)::ConvexObjective{Th,Tg} where {Th<:Union{Matrix{Float},SpMatrix{Float}},Tg<:Union{Vector{Float},SpVector{Float}}}

    # collect info
    numVars = get_numVariables(objective)

    L = Tg(objective.L)
    c = objective.c
    if Th == Matrix{Float}
        H = zeros(numVars,numVars)
    else
        H = spzeros(numVars,numVars)
    end

    # define value function, gradient and hessian
    newEvalVal = x::Vector{Float}->(L'*x + c)::Float
    newEvalGrd = x::Vector{Float}->L::Tg
    newEvalHes = x::Vector{Float}->H::SpMatrix{Float}
    return ConvexObjective{Th,Tg}(newEvalVal,newEvalGrd,newEvalHes,sparsityVector(objective.L),spzeros(Bool,numVars,numVars))
end

function SparseArrays.sparse(objective::ConvexObjective)::ConvexObjective{SpMatrix{Float},SpVector{Float}}
    return ConvexObjective{SpMatrix{Float},SpVector{Float}}(objective)
end

# inspect functions  (Fundamental. These are used in Branch and Bound)
function get_numVariables(objective::ConvexObjective)::Int
    return length(objective.grdSparsity)
end

function get_dependency(objective::ConvexObjective)::Vector{Int}
    return findnz(objective.grdSparsity)[1]
end

function get_gradientSparsity(objective::ConvexObjective)::Vector{Int}
    return findnz(objective.grdSparsity)[1]
end

function get_hessianSparsity(objective::ConvexObjective)::Tuple{Vector{Int},Vector{Int}}
    return findnz(objective.hesSparsity)[1:2]
end

function get_gradientType(constraintSet::ConvexObjective{<:AbstractMatrix,Tg})::Type where Tg<:Union{Vector{Float},SpVector{Float}}
    return Tg
end

function get_hessianType(constraintSet::ConvexObjective{Th})::Type where Th<:Union{Matrix{Float},SpMatrix{Float}}
    return Th
end

# evaluate the objective function in the given point
function evaluate(objective::ConvexObjective,point::VirtualVector{Float})::Float
    @assert length(point) == get_numVariables(objective)
    return objective.evalVal(point)
end

# evaluate the gradient of the objective function in the given point
function evaluate_gradient(objective::ConvexObjective{Th,Tg},point::VirtualVector{Float})::Tg where {Th<:Union{Matrix{Float},SpMatrix{Float}},Tg<:Union{Vector{Float},SpVector{Float}}}
    @assert length(point) == get_numVariables(objective)
    return objective.evalGrd(point)
end

# evaluate the hessian of the objective function in the given point
function evaluate_hessian(objective::ConvexObjective{Th,Tg},point::VirtualVector{Float})::Th where {Th<:Union{Matrix{Float},SpMatrix{Float}},Tg<:Union{Vector{Float},SpVector{Float}}}
    @assert length(point) == get_numVariables(objective)
    return objective.evalHes(point)
end

# update functions (Not fundamental These are used only for problem update)
function insert_variables!(objective::ConvexObjective{Th,Tg},numVariables::Int,insertionPoint::Int)::Nothing where {Th<:Union{Matrix{Float},SpMatrix{Float}},Tg<:Union{Vector{Float},SpVector{Float}}}
    @assert numVariables >= 0
    @assert 1<=insertionPoint<=get_numVariables(objective)+1

    # collect info
    indices = append(collect(1:insertionPoint-1),collect(insertionPoint+numVariables:oldNumVars))

    # define the new value function (wrapping the old one)
    oldEvalVal = objective.evalVal
    objective.evalVal = x::Vector{Float}->oldEvalVal(x[indices])::Float

    # define new gradient
    if Tg == Vector{Float}
        oldEvalGrd = objective.evalGrd
        function newEvalGrdD(x::Vector{Float})::Vector{Float}
            gradient_ = zeros(length(x))
            gradient_[indices] = oldEvalGrd(x[indices])
            return gradient_
        end
        objective.evalGrd = newEvalGrdD
    else
        oldEvalGrd = objective.evalGrd
        function newEvalGrdS(x::Vector{Float})::SpVector{Float}
            gradient_ = spzeros(length(x))
            gradient_[indices] = oldEvalGrd(x[indices])
            return gradient_
        end
        objective.evalGrd = newEvalGrdS
    end


    # define new hessian
    if Th == Matrix{Float}
        oldEvalHes = objective.evalHes
        function newEvalHesD(x::Vector{Float})::Matrix{Float}
            hessian_ = zeros(newNumVars,newNumVars)
            hessian_[indices,indices] = oldEvalHes(x[indices])
            return hessian_
        end
        objective.evalHes = newEvalHesD
    else
        oldEvalHes = objective.evalHes
        function newEvalHesS(x::Vector{Float})::SpMatrix{Float}
            hessian_ = spzeros(newNumVars,newNumVars)
            hessian_[indices,indices] = oldEvalHes(x[indices])
            return hessian_
        end
        objective.evalHes = newEvalHesS
    end


    # update gradient sparsity
    grdSparsity_ = sparsevec(Int[],Bool[],newNumVars)
    grdSparsity_[indices] = objective.grdSparsity
    objective.grdSparsity = grdSparsity_

    # update hessian sparsity
    hesSparsity_ = sparse(Int[],Int[],Bool[],newNumVars,newNumVars)
    hesSparsity_[indices,indices] = objective.hesSparsity
    objective.hesSparsity = hesSparsity_

    return
end

function append_variables!(objective::ConvexObjective,numVariables::Int)::Nothing
    insert_variables!(objective,numVariables,get_numVariables(objective)+1)
    return
end


function remove_variables!(objective::ConvexObjective{Th,Tg},indices::Vector{Int})::Nothing where {Th<:Union{Matrix{Float},SpMatrix{Float}},Tg<:Union{Vector{Float},SpVector{Float}}}

    # check correctness of input
    if !iszero(objective.grdSparsity([indices]))
        error("The objective function depends on the variables to be removed: impossible to remove them")
    end

    # collect info
    oldNumVars = get_numVariables(objective)
    toKeep = filter(x->!(x in indices),collect(1:oldNumVars))


    # define the new value function (wrapping the old one)
    oldEvalVal = objective.evalVal
    x_ = Vector{Float}(undef,oldNumVars)
    function newEvalVal(x::Vector{Float})::Float
        x_[toKeep] = x
        return oldEvalVal(x_)
    end
    objective.evalVal = newEvalVal

    # update gradient
    if Tg == Vector{Float}
        oldEvalGrd = objective.evalGrd
        x_ = Vector{Float}(undef,oldNumVars)
        function newEvalGrdD(x::Vector{Float})::Vector{Float}
            x_[toKeep] = x
            return oldEvalGrd(x_)[toKeep]
        end
        objective.evalGrd = newEvalGrdD
    else
        oldEvalGrd = objective.evalGrd
        x_ = Vector{Float}(undef,oldNumVars)
        function newEvalGrdS(x::Vector{Float})::SpVector{Float}
            x_[toKeep] = x
            return oldEvalGrd(x_)[toKeep]
        end
        objective.evalGrd = newEvalGrdS
    end

    # update hessian
    if Th == Matrix{Float}
        oldEvalHes = objective.evalHes
        x_ = Vector{Float}(undef,oldNumVars)
        function newEvalHesD(x::Vector{Float})::Matrix{Float}
            x_[toKeep] = x
            return oldEvalHes(x_)[toKeep,toKeep]
        end
        objective.evalHes = newEvalHesD
    else
        oldEvalHes = objective.evalHes
        x_ = Vector{Float}(undef,oldNumVars)
        function newEvalHesS(x::Vector{Float})::SpMatrix{Float}
            x_[toKeep] = x
            return oldEvalHes(x_)[toKeep,toKeep]
        end
        objective.evalHes = newEvalHesS
    end

    # update sparsities
    objective.grdSparsity = objective.grdSparsity[toKeep]
    objective.hesSparsity = objective.hesSparsity[toKeep,toKeep]

    return
end



function fix_variables!(objective::ConvexObjective{Th,Tg},indices::Vector{Int},values::Vector{Float};removeFixedVariables::Bool=false)::Nothing where {Th<:Union{Matrix{Float},SpMatrix{Float}},Tg<:Union{Vector{Float},SpVector{Float}}}

    # collect info
    oldNumVars = get_numVariables(objective)
    toKeep = filter(x->!(x in indices),collect(1:oldNumVars))


    # define the new value function (wrapping the old one)
    oldEvalVal = objective.evalVal
    x_ = Vector{Float}(undef,oldNumVars); x_[indices] = values
    if removeFixedVariables
        function newEvalVal1(x::Vector{Float})::Float
            x_[toKeep] = x
            return oldEvalVal(x_)
        end
        objective.evalVal = newEvalVal1
    else
        function newEvalVal2(x::Vector{Float})::Float
            x_[toKeep] = x[toKeep]
            return oldEvalVal(x_)
        end
        objective.evalVal = newEvalVal2
    end

    # update gradient
    if Tg == Vector{Float}
        oldEvalGrd = objective.evalGrd
        x_ = Vector{Float}(undef,oldNumVars); x_[indices] = values
        if removeFixedVariables
            function newEvalGrdD1(x::Vector{Float})::Vector{Float}
                x_[toKeep] = x
                return oldEvalGrd(x_)[toKeep]
            end
            objective.evalGrd = newEvalGrdD1
        else
            function newEvalGrdD2(x::Vector{Float})::Vector{Float}
                x_[toKeep] = x[toKeep]
                gradient_ = oldEvalGrd(x_)
                gradient_[indices] = zeros(length(indices))
                return gradient_
            end
            objective.evalGrd = newEvalGrdD2
        end
    else
        oldEvalGrd = objective.evalGrd
        x_ = Vector{Float}(oldNumVars); x_[indices] = values
        if removeFixedVariables
            function newEvalGrdS1(x::Vector{Float})::SpVector{Float}
                x_[toKeep] = x
                return oldEvalGrd(x_)[toKeep]
            end
            objective.evalGrd = newEvalGrdS1
        else
            function newEvalGrdS2(x::Vector{Float})::SpVector{Float}
                x_[toKeep] = x[toKeep]
                gradient_ = oldEvalGrd(x_)
                gradient_[indices] = spzeros(length(indices))
                return gradient_
            end
            objective.evalGrd = newEvalGrdS2
        end
    end

    # update hessian
    if Th == Matrix{Float}
        oldEvalHes = objective.evalHes
        x_ = Vector{Float}(undef,oldNumVars); x_[indices] = values
        if removeFixedVariables
            function newEvalHesD1(x::Vector{Float})::Matrix{Float}
                x_[toKeep] = x
                return oldEvalHes(x_)[toKeep,toKeep]
            end
            objective.evalHes = newEvalHesD1
        else
            function newEvalHesD2(x::Vector{Float})::Matrix{Float}
                x_[toKeep] = x[toKeep]
                hessian_ = oldEvalHes(x_)
                hessian_[:,indices] = zeros(size(hessian_,1),length(indices))
                hessian_[indices,:] = zeros(length(indices),size(hessian_,2))
                return hessian_
            end
            objective.evalHes = newEvalHesD2
        end
    else
        oldEvalHes = objective.evalHes
        x_ = Vector{Float}(undef,oldNumVars); x_[indices] = values
        if removeFixedVariables
            function newEvalHesS1(x::Vector{Float})::SpMatrix{Float}
                x_[toKeep] = x
                return oldEvalHes(x_)[toKeep,toKeep]
            end
            objective.evalHes = newEvalHesS1
        else
            function newEvalHesS2(x::Vector{Float})::SpMatrix{Float}
                x_[toKeep] = x
                hessian_[:,indices] = spzeros(size(hessian_,1),length(indices))
                hessian_[indices,:] = spzeros(length(indices),size(hessian_,2))
                return hessian_
            end
            objective.evalHes = newEvalHesS2
        end
    end

    # update sparsities
    if removeFixedVariables
        objective.grdSparsity = objective.grdSparsity[toKeep]
        objective.hesSparsity = objective.hesSparsity[toKeep,toKeep]
    else
        objective.grdSparsity[indices] = spzeros(Bool,length(indices))
        objective.hesSparsity[:,indices] = spzeros(Bool,size(objective.hesSparsity,1),length(indices))
        objective.hesSparsity[indices,:] = spzeros(Bool,length(indices),size(objective.hesSparsity,2))
    end

    return
end

import Base.+
function +(objective1::ConvexObjective{Th,Tg},objective2::ConvexObjective{Th,Tg})::ConvexObjective{Th,Tg} where {Th<:Union{Matrix{Float},SpMatrix{Float}},Tg<:Union{Vector{Float},SpVector{Float}}}

    # build the new functions
    (V1,V2) = (objective1.evalVal,objective2.evalVal)
    function newEvalVal(x::Vector{Float})::Float
        return V1(x) + V2(x)
    end
    evalVal = newEvalVal

    if Tg == Vector{Float}
        (G1,G2) = (objective1.evalGrd,objective2.evalGrd)
        function newEvalGrdD(x::Vector{Float})::Vector{Float}
            return G1(x) + G2(x)
        end
        evalGrd = newEvalGrdD
    else
        (G1,G2) = (objective1.evalGrd,objective2.evalGrd)
        function newEvalGrdS(x::Vector{Float})::SpVector{Float}
            return G1(x) + G2(x)
        end
        evalGrd = newEvalGrdS
    end


    if Th == Matrix{Float}
        (H1,H2) = (objective1.evalHes,objective2.evalHes)
        function newEvalHesD(x::Vector{Float})::Matrix{Float}
            return H1(x) + H2(x)
        end
        evalHes = newEvalHesD
    else
        (H1,H2) = (objective1.evalHes,objective2.evalHes)
        function newEvalHesS(x::Vector{Float})::SpMatrix{Float}
            return H1(x) + H2(x)
        end
        evalHes = newEvalHesS
    end

    grdSparsity = dropzeros(objective1.grdSparsity .| objective2.grdSparsity)
    hesSparsity = dropzeros(objective1.hesSparsity .| objective2.hesSparsity)

    return ConvexObjective{Th,Tg}(evalVal,evalHes,evalGrd,grdSparsity,hesSparsity)
end


function add!(objective1::ConvexObjective{Th,Tg},objective2::ConvexObjective{Th,Tg})::Nothing where {Th<:Union{Matrix{Float},SpMatrix{Float}},Tg<:Union{Vector{Float},SpVector{Float}}}

    # build the new functions
    (V1,V2) = (objective1.evalVal,objective2.evalVal)
    function newEvalVal(x::Vector{Float})::Float
        return V1(x) + V2(x)
    end
    objective1.evalVal = newEvalVal

    if Tg == Vector{Float}
        (G1,G2) = (objective1.evalGrd,objective2.evalGrd)
        function newEvalGrdD(x::Vector{Float})::Vector{Float}
            return G1(x) + G2(x)
        end
        objective1.evalGrd = newEvalGrdD
    else
        (G1,G2) = (objective1.evalGrd,objective2.evalGrd)
        function newEvalGrdS(x::Vector{Float})::SpVector{Float}
            return G1(x) + G2(x)
        end
        objective1.evalGrd = newEvalGrdS
    end


    if Th == Matrix{Float}
        (H1,H2) = (objective1.evalHes,objective2.evalHes)
        function newEvalHesD(x::Vector{Float})::Matrix{Float}
            return H1(x) + H2(x)
        end
        objective1.evalHes = newEvalHesD
    else
        (H1,H2) = (objective1.evalHes,objective2.evalHes)
        function newEvalHesS(x::Vector{Float})::SpMatrix{Float}
            return H1(x) + H2(x)
        end
        objective1.evalHes = newEvalHesS
    end

    objective1.grdSparsity = dropzeros(objective1.grdSparsity .| objective2.grdSparsity)
    objective1.hesSparsity = dropzeros(objective1.hesSparsity .| objective2.hesSparsity)

    return
end






for type in ConvexObjectiveLeafTypes @assert precompile(is_linear,(type,)) end
for type in ConvexObjectiveLeafTypes @assert precompile(is_quadratic,(type,)) end
for type in ConvexObjectiveLeafTypes @assert precompile(is_convex,(type,)) end
for type in ConvexObjectiveLeafTypes @assert precompile(copy,(type,)) end
for type in ConvexObjectiveLeafTypes @assert precompile(deepcopy,(type,)) end
for type1 in ConvexObjectiveLeafTypes
    for type2 in ConvexObjectiveLeafTypes @assert precompile(type1,(type2,)) end
    for type2 in QuadraticObjectiveLeafTypes @assert precompile(type1,(type2,)) end
    for type2 in LinearObjectiveLeafTypes @assert precompile(type1,(type2,)) end
end
for type in ConvexObjectiveLeafTypes @assert precompile(sparse,(type,)) end
for type in ConvexObjectiveLeafTypes @assert precompile(get_numVariables,(type,)) end
for type in ConvexObjectiveLeafTypes @assert precompile(get_dependency,(type,)) end
for type in ConvexObjectiveLeafTypes @assert precompile(get_gradientSparsity,(type,)) end
for type in ConvexObjectiveLeafTypes @assert precompile(get_hessianSparsity,(type,)) end
for type in ConvexObjectiveLeafTypes @assert precompile(get_gradientType,(type,)) end
for type in ConvexObjectiveLeafTypes @assert precompile(get_hessianType,(type,)) end
for type in ConvexObjectiveLeafTypes @assert precompile(firstOrderTaylor,(type,Vector{Float})) end
for type in ConvexObjectiveLeafTypes @assert precompile(evaluate,(type,Vector{Float})) end
for type in ConvexObjectiveLeafTypes @assert precompile(evaluate_gradient,(type,Vector{Float})) end
for type in ConvexObjectiveLeafTypes @assert precompile(evaluate_hessian,(type,Vector{Float})) end
