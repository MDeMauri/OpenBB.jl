# @Author: Massimo De Mauri <massimo>
# @Date:   2019-05-23T11:28:46+02:00
# @Email:  massimo.demauri@gmail.com
# @Filename: BBstatus.jl
# @Last modified by:   massimo
# @Last modified time: 2019-08-12T22:04:25+02:00
# @License: LGPL-3.0
# @Copyright: {{copyright}}


mutable struct BBstatus
    objLoB::Float64
    objUpB::Float64
    absoluteGap::Float64
    relativeGap::Float64
    totalTime::Float64
    waitingTime::Float64
    numSolutions::Int
    numExploredNodes::Int64
    reliable::Bool
    cutoffActive::Bool
    description::String
end

function BBstatus(; objLoB::Float64=-Inf,objUpB::Float64=Inf,
                    absoluteGap::Float64= Inf,relativeGap::Float64=Inf,
                    totalTime::Float64=0.0,waitingTime::Float64=0.0,
                    numSolutions::Int=0,numExploredNodes::Int=0,
                    reliable::Bool=true,
                    cutoffActive::Bool=false,
                    description::String="new")::BBstatus

    return BBstatus(objLoB,objUpB,absoluteGap,relativeGap,
                    totalTime,waitingTime,
                    numSolutions,numExploredNodes,
                    reliable,cutoffActive,description)
end


function BBstatus(status::BBstatus)::BBstatus

    return BBstatus(status.objLoB,status.objUpB,
                    status.absoluteGap,status.relativeGap,
                    status.totalTime,status.waitingTime,
                    status.numSolutions,status.numExploredNodes,
                    status.reliable,status.cutoffActive,status.description)
end
