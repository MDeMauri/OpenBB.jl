# @Author: Massimo De Mauri <massimo>
# @Date:   2019-02-26T15:03:06+01:00
# @Email:  massimo.demauri@gmail.com
# @Filename: ObjectiveFunctions.jl
# @Last modified by:   massimo
# @Last modified time: 2021-01-16T16:06:25+01:00
# @License: LGPL-3.0
# @Copyright: {{copyright}}

# abstract and null types
abstract type ObjectiveFunction end
struct NullObjective <: ObjectiveFunction end

## Generic Functions
import Base.+
function +(objFun1::T1,objFun2::T2)::T1 where {T1<:ObjectiveFunction,T2<:ObjectiveFunction}
    return objFun1 + T1(objFun2)
end

function add!(objFun1::T1,objFun2::T2)::Nothing where {T1<:ObjectiveFunction,T2<:ObjectiveFunction}
    add!(objFun1,T1(objFun2))
    return
end

# linear approximation (relaxation if convex) of the objective in the given point
function firstOrderTaylor(objFun::ObjectiveFunction,point::VirtualVector{Float})::LinearObjective
    L = evaluate_gradient(objFun,point)
    c = evaluate(objFun,point) - L'*point
    return LinearObjective(L=L,c=c)
end

######################## Type Definitions ##############################
mutable struct LinearObjective{T<:Union{Vector{Float},SpVector{Float}}} <: ObjectiveFunction
    L::T
    c::Float
end
LinearObjectiveLeafTypes = [LinearObjective{Vector{Float}},LinearObjective{SpVector{Float}}]


############
mutable struct QuadraticObjective{TQ<:Union{Matrix{Float},SpMatrix{Float}},TL<:Union{Vector{Float},SpVector{Float}}} <: ObjectiveFunction
    Q::TQ
    L::TL
    c::Float
    pInvQ::Union{Matrix{Float},Nothing}
end
QuadraticObjectiveLeafTypes = [QuadraticObjective{Matrix{Float},Vector{Float}},
                               QuadraticObjective{Matrix{Float},SpVector{Float}},
                               QuadraticObjective{SpMatrix{Float},Vector{Float}},
                               QuadraticObjective{SpMatrix{Float},SpVector{Float}}]


#############
mutable struct ConvexObjective{Th<:Union{Matrix{Float},SpMatrix{Float}},Tg<:Union{Vector{Float},SpVector{Float}}} <: ObjectiveFunction
    evalVal::Function
    evalGrd::Function
    evalHes::Function
    grdSparsity::SpVector{Bool}
    hesSparsity::SpMatrix{Bool}

    # non-standard constructor: Implements a probing step to ensure the correctness of the function definitions
    function ConvexObjective{Th,Tg}(evalVal::Function,evalGrd::Function,evalHes::Function,
                                       grdSparsity::SpVector{Bool},hesSparsity::SpMatrix{Bool};probingPoint::Vector{Float}=Float[]
                                       )::ConvexObjective{Th,Tg} where {Th<:Union{Matrix{Float},SpMatrix{Float}},Tg<:Union{Vector{Float},SpVector{Float}}}

        # detect num variables
        numVars = length(grdSparsity)
        # generate a probing point if none is given
        if isempty(probingPoint)
            probingPoint = rand(numVars)
        end
        # probing (to avoid future problems with the input output types of the functions defining the objective)
        try
            evalVal(probingPoint)::Float
            evalGrd(probingPoint)::Tg
            evalHes(probingPoint)::Th
        catch err
            if err isa MethodError
                @info err
                throw(ArgumentError("The Input types of the functions defining the objective are not all correct"))
            elseif err isa TypeError
                @info err
                throw(ArgumentError("The output types of the functions defining the objective are not all correct"))
            else
                @info err
                @warn "An unknown error while checking the input-output types the functions defining the objective.
                       This may mean that the used probing point hit an undefined point for one of the functions, e.g. log(-1),
                       or some other unimportant problem arised. However, check input-output dimensions and types to be sure."
            end
        end

        # wrap the functions for allowing the use of views (splicing of array without copying)
        x_ = Vector{Float}(undef,numVars)
        function wrapEvalVal(x::VirtualVector{Float})::Float
            if x isa SubArray x_ .= x; return evalVal(x_) else return evalVal(x) end
        end
        function wrapEvalGrd(x::VirtualVector{Float})::Tg
            if x isa SubArray x_ .= x; return evalGrd(x_) else return evalGrd(x) end
        end
        function wrapEvalHes(x::VirtualVector{Float})::Th
            if x isa SubArray x_ .= x; return evalHes(x_) else return evalHes(x) end
        end

        # actually populate the struct
        return new{Th,Tg}(wrapEvalVal,wrapEvalGrd,wrapEvalHes,grdSparsity,hesSparsity)
    end
end
ConvexObjectiveLeafTypes = [ConvexObjective{Matrix{Float},Vector{Float}},
                            ConvexObjective{Matrix{Float},SpVector{Float}},
                            ConvexObjective{SpMatrix{Float},Vector{Float}},
                            ConvexObjective{SpMatrix{Float},SpVector{Float}}]
