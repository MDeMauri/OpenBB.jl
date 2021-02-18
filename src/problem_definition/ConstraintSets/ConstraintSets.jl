# @Author: Massimo De Mauri <massimo>
# @Date:   2019-02-26T15:03:06+01:00
# @Email:  massimo.demauri@gmail.com
# @Filename: ConstraintSets.jl
# @Last modified by:   massimo
# @Last modified time: 2021-02-12T14:44:25+01:00
# @License: LGPL-3.0
# @Copyright: {{copyright}}

# Note: here we use non-standard types defined in OpenBB,jl

# abstract and null types
abstract type ConstraintSet end
struct NullConstraintSet <: ConstraintSet end

## common operations
function Base.insert!(constraintSet1::T1,constraintSet2::T2,insertionPoint::Int)::Nothing where {T1<:ConstraintSet,T2<:ConstraintSet}
    insert!(constraintSet1,T1(constraintSet2))
end

function Base.append!(constraintSet1::T1,constraintSet2::T2)::Nothing where {T1<:ConstraintSet,T2<:ConstraintSet}
    insert!(constraintSet1,T1(constraintSet2),get_size(constraintSet1)+1)
end

function concat(constraintSet1::T1,constraintSet2::T2)::T1 where {T1<:ConstraintSet,T2<:ConstraintSet}
    constraintSet_ = deepcopy(constraintSet1)
    insert!(constraintSet_,T1(constraintSet2),get_size(constraintSet1)+1)
    return constraintSet_
end

function Base.getindex(constraintSet::T,list::Vector{I})::T where {T<:ConstraintSet,I<:Integer}
    toRemove = [k for k in 1:get_size(constraintSet) if !(k in list)]
    newConstraintSet = deepcopy(constraintSet)
    remove_constraints!(newConstraintSet,toRemove)
    return newConstraintSet
end

function Base.getindex(constraintSet::T,range::UnitRange)::T where {T<:ConstraintSet}
    return getindex(constraintSet,collect(range))
end

function Base.getindex(constraintSet::T,index::I)::T where {T<:ConstraintSet,I<:Integer}
    return getindex(constraintSet,[index])
end


function sort_oc!(numVarsPerStep::Int,constraintSet::T)::Nothing where T <: ConstraintSet
    # get the first non-zero element (row-wise) of each constraint
    timeSteps = div.(get_firstNZs(workspace.problem.cnsSet,2).+(numVarsPerStep-1),numVarsPerStep)
    # sort the constraints
    if !issorted(timeSteps)
        permute!(constraintSet,sortperm(timeSteps,alg=MergeSort)) # merge sort for stability (not the fastest)
    end
    return
end

## Optimal Control Specific functions

# identify the time steps the constraints belong to assuming a
# single/multiple shooting structure of the problem (also collocation?)
function get_timeSteps_oc(cnsSet::ConstraintSet,optimalControlInfo::Tuple{Int,Int,Int})::Vector{Int}
    @assert all(optimalControlInfo .>=0)
    numVarsPerStep = sum(optimalControlInfo)
    return div.((numVarsPerStep-1) .+ get_firstNZs(cnsSet,2),numVarsPerStep)
end


function issorted_oc(cnsSet::ConstraintSet,optimalControlInfo::Tuple{Int,Int,Int})::Bool
    return issorted(get_timeSteps_oc(cnsSet,optimalControlInfo))
end

function sortperm_oc(cnsSet::ConstraintSet,optimalControlInfo::Tuple{Int,Int,Int})::Vector{Int}
    if !issorted_oc(cnsSet,optimalControlInfo)
        return sortperm(get_timeSteps_oc(cnsSet,optimalControlInfo),alg=MergeSort)
    else
        return collect(1:get_size(cnsSet))
    end
end




############################## Type Definitions ########################3
mutable struct LinearConstraintSet{T} <: ConstraintSet where T<:Union{Matrix{Float},SpMatrix{Float}}
    A::T
    loBs::Vector{Float}
    upBs::Vector{Float}
end
LinearConstraintSetLeafTypes = [LinearConstraintSet{Matrix{Float}},
                                LinearConstraintSet{SpMatrix{Float}}]

#################
mutable struct ConvexConstraintSet{Th,Tj} <: ConstraintSet where {Th<:Union{Matrix{Float},SpMatrix{Float}},Tj<:Union{Matrix{Float},SpMatrix{Float}}}
    evalVal::Function
    evalJcb::Function
    evalHes::Function
    jcbSparsity::SpMatrix{Bool}
    hesSparsity::Vector{SpMatrix{Bool}}
    loBs::Vector{Float}
    upBs::Vector{Float}

    # non-standard constructor: Implements a probing step to ensure the correctness of the input-output types of the provided function definitions
    function ConvexConstraintSet{Th,Tj}(evalVal::Function,evalJcb::Function,evalHes::Function,
                                           jcbSparsity::SpMatrix{Bool},hesSparsity::Vector{SpMatrix{Bool}},
                                           loBs::Vector{Float},upBs::Vector{Float};probingPoint::Vector{Float}=Float[]
                                          )::ConvexConstraintSet{Th,Tj} where {Th<:Union{Matrix{Float},SpMatrix{Float}},Tj<:Union{Matrix{Float},SpMatrix{Float}}}

        # detect num variables
        numVars = size(jcbSparsity,2)
        # generate a probing point if none is provided
        if isempty(probingPoint)
            probingPoint=rand(numVars)
        end
        # probing (to avoid future problems with the input output types of the functions defining the constraint set)
        try
          evalVal(probingPoint)::Vector{Float}
          evalJcb(probingPoint)::Tj
          evalHes(probingPoint)::Vector{Th}
        catch err
          if typeof(err) isa MethodError
              error("The Input types of the functions defining the constraint set are not all correct")
          elseif typeof(err) isa TypeError
              error("The output types of the functions defining the constraint set are not all correct")
          else
              @info err
              @warn "An unknown error while checking the input-output types the function.
                     This may mean that the used probing point hit an undefined point for the function, e.g. log(-1),
                     or some other unimportant problem arised. However, check dimensions and types to be sure."
          end
        end

        # wrap the functions for allowing the use of views (splicing of array without copying)
        x_ = Vector{Float}(undef,numVars)
        function wrapEvalVal(x::VirtualVector{Float})::Vector{Float}
            if x isa SubArray x_ .= x; return evalVal(x_) else return evalVal(x) end
        end
        function wrapEvalJcb(x::VirtualVector{Float})::Tj
            if x isa SubArray x_ .= x; return evalJcb(x_) else return evalJcb(x) end
        end
        function wrapEvalHes(x::VirtualVector{Float})::Vector{Th}
            if x isa SubArray x_ .= x; return evalHes(x_) else return evalHes(x) end
        end

        return new{Th,Tj}(wrapEvalVal,wrapEvalJcb,wrapEvalHes,jcbSparsity,hesSparsity,loBs,upBs)
    end
end
ConvexConstraintSetLeafTypes = [ConvexConstraintSet{Matrix{Float},Matrix{Float}},
                                ConvexConstraintSet{Matrix{Float},SpMatrix{Float}},
                                ConvexConstraintSet{SpMatrix{Float},Matrix{Float}},
                                ConvexConstraintSet{SpMatrix{Float},SpMatrix{Float}}]
