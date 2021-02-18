# @Author: Massimo De Mauri <massimo>
# @Date:   2019-05-23T11:28:46+02:00
# @Email:  massimo.demauri@gmail.com
# @Filename: BBstatus.jl
# @Last modified by:   massimo
# @Last modified time: 2020-12-09T13:37:36+01:00
# @License: LGPL-3.0
# @Copyright: {{copyright}}


mutable struct BBstatus
    objLoB::Float
    objUpB::Float
    absoluteGap::Float
    relativeGap::Float
    totalTime::Float
    waitingTime::Float
    numSolutions::Int
    numExploredNodes::Int64
    reliable::Bool
    cutoffActive::Bool
    blackListActive::Bool
    somethingWrong::Bool
    description::String
end

function BBstatus(; objLoB::Float=-Inf,objUpB::Float=Inf,
                    absoluteGap::Float= Inf,relativeGap::Float=Inf,
                    totalTime::Float=0.0,waitingTime::Float=0.0,
                    numSolutions::Int=0,numExploredNodes::Int=0,
                    reliable::Bool=true,
                    cutoffActive::Bool=false,
                    blackListActive::Bool=false,
                    somethingWrong::Bool=false,
                    description::String="new")::BBstatus

    return BBstatus(objLoB,objUpB,absoluteGap,relativeGap,
                    totalTime,waitingTime,
                    numSolutions,numExploredNodes,
                    reliable,cutoffActive,blackListActive,
                    somethingWrong,description)
end



function Base.copy(status::BBstatus)::BBstatus

    return BBstatus(status.objLoB,status.objUpB,
                    status.absoluteGap,status.relativeGap,
                    status.totalTime,status.waitingTime,
                    status.numSolutions,status.numExploredNodes,
                    status.reliable,status.cutoffActive,status.blackListActive,
                    status.somethingWrong,status.description)
end

function Base.deepcopy(status::BBstatus)::BBstatus

    return BBstatus(copy(status.objLoB),copy(status.objUpB),
                    status.absoluteGap,status.relativeGap,
                    status.totalTime,status.waitingTime,
                    status.numSolutions,status.numExploredNodes,
                    status.reliable,status.cutoffActive,status.blackListActive,
                    status.somethingWrong,status.description)
end
