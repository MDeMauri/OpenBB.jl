# @Author: Massimo De Mauri <massimo>
# @Date:   2020-11-23T11:47:15+01:00
# @Email:  massimo.demauri@gmail.com
# @Filename: HBBstatus.jl
# @Last modified by:   massimo
# @Last modified time: 2020-11-26T18:00:05+01:00
# @License: LGPL-3.0
# @Copyright: {{copyright}}

mutable struct HBBstatus
    objLoB::Float
    objUpB::Float
    absoluteGap::Float
    relativeGap::Float
    totalTime::Float
    nlpTime::Float
    mipTime::Float
    rlxTime::Float
    numIncumbents::Int
    numIterations::Int
    reliable::Bool
    cutoffActive::Bool
    description::String
end


function HBBstatus(;objLoB::Float=-Inf,objUpB::Float=Inf,
                    absoluteGap::Float=Inf,relativeGap::Float=Inf,
                    totalTime::Float=0.0,nlpTime::Float=0.0,
                    mipTime::Float=0.0,rlxTime::Float=0.0,
                    numIncumbents::Int=0,numIterations::Int=0,
                    reliable::Bool=true,
                    cutoffActive::Bool=false,
                    description::String="new")::HBBstatus

    return HBBstatus(objLoB,objUpB,absoluteGap,relativeGap,
                     totalTime,nlpTime,mipTime,rlxTime,
                     numIncumbents,numIterations,
                     reliable,cutoffActive,description)
end


function Base.copy(status::HBBstatus)::HBBstatus

    return HBBstatus(status.objLoB,status.objUpB,
                     status.absoluteGap,status.relativeGap,
                     status.totalTime,status.nlpTime,status.mipTime,status.rlxTime,
                     status.numIncumbents,status.numIterations,
                     status.reliable,status.cutoffActive,status.description)
end

function Base.deepcopy(status::HBBstatus)::HBBstatus

    return HBBstatus(copy(status.objLoB),copy(status.objUpB),
                    status.absoluteGap,status.relativeGap,
                    status.totalTime,status.nlpTime,status.mipTime,status.rlxTime,
                    status.numIncumbents,status.numIterations,
                    status.reliable,status.cutoffActive,status.description)
end
