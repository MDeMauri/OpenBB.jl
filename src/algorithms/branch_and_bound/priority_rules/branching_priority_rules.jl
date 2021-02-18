# @Author: Massimo De Mauri <massimodemauri>
# @Date:   2019-02-06T16:19:41+01:00
# @Email:  massimo.demauri@gmail.com
# @Project: OpenBB
# @Filename: branching_priority_functions.jl
# @Last modified by:   massimo
# @Last modified time: 2021-01-08T19:28:07+01:00
# @License: LGPL-3.0
# @Copyright: {{copyright}}


# wrapper for branching priority rules
function branching_priority_rule(functionTuple::Tuple{Symbol,Vararg},dscValues::Vector{Float},pseudoCosts::Tuple{Matrix{Float},Matrix{Int}},primalTolerance::Float)::Tuple{Int64,Float}
    return OpenBB.eval(functionTuple[1])(dscValues,pseudoCosts,primalTolerance,functionTuple[2:end]...)
end

# actual priority rules
function pseudoIncrements_mean(dscValues::Vector{Float},pseudoCosts::Tuple{Matrix{Float},Matrix{Int}},primalTolerance::Float)::Tuple{Int64,Float}
    return pseudoIncrements_mean(dscValues,pseudoCosts,primalTolerance,1,1)
end

function pseudoIncrements_mean(dscValues::Vector{Float},pseudoCosts::Tuple{Matrix{Float},Matrix{Int}},primalTolerance::Float,
                          ratio::Rational)::Tuple{Int64,Float}
    return pseudoIncrements_mean(dscValues,pseudoCosts,primalTolerance,numerator(ratio),denominator(ratio)-numerator(ratio))
end

function pseudoIncrements_mean(dscValues::Vector{Float},pseudoCosts::Tuple{Matrix{Float},Matrix{Int}},primalTolerance::Float,
                               gainOfMin::Int64,gainOfMax::Int64)::Tuple{Int64,Float}

    bestIndex = 0
    bestScore = 0.0
    for k in 1: length(dscValues)
        deltaMinus = threshold(dscValues[k] - floor(dscValues[k]),primalTolerance)*pseudoCosts[1][k,1]
        deltaPlus  = threshold(ceil(dscValues[k]) - dscValues[k],primalTolerance)*pseudoCosts[1][k,2]
        if deltaMinus <= deltaPlus
            score_ = deltaMinus*gainOfMin + deltaPlus*gainOfMax
        else
            score_ = deltaPlus*gainOfMin + deltaMinus*gainOfMax
        end

        if score_ > bestScore
            bestIndex = k
            bestScore = score_
        end
    end
    return bestIndex, bestScore
end

function pseudoIncrements_geomean(dscValues::Vector{Float},pseudoCosts::Tuple{Matrix{Float},Matrix{Int}},primalTolerance::Float)::Tuple{Int64,Float}
    return pseudoIncrements_geomean(dscValues,pseudoCosts,primalTolerance,1,1)
end

function pseudoIncrements_geomean(dscValues::Vector{Float},pseudoCosts::Tuple{Matrix{Float},Matrix{Int}},primalTolerance::Float,
                                      ratio::Rational)::Tuple{Int64,Float}
    return pseudoIncrements_geomean(dscValues,pseudoCosts,primalTolerance,numerator(ratio),denominator(ratio)-numerator(ratio))
end

function pseudoIncrements_geomean(dscValues::Vector{Float},pseudoCosts::Tuple{Matrix{Float},Matrix{Int}},primalTolerance::Float,
                             gainOfMin::Int,gainOfMax::Int)::Tuple{Int64,Float}

    score_s = zeros(length(dscValues))
    bestIndex = 0
    bestScore = 0.0
    for k in 1:length(dscValues)
        deltaMinus = threshold(dscValues[k] - floor(dscValues[k]),primalTolerance)*pseudoCosts[1][k,1]
        deltaPlus  = threshold(ceil(dscValues[k]) - dscValues[k],primalTolerance)*pseudoCosts[1][k,2]

        if deltaMinus <= deltaPlus
            score_ = deltaMinus^gainOfMin*deltaPlus^gainOfMax
        else
            score_ = deltaPlus^gainOfMin*deltaMinus^gainOfMax
        end

        if score_ > bestScore
            bestIndex = k
            bestScore = score_
        end
    end
    return bestIndex, bestScore
end



function absolute_fractionality(dscValues::Vector{Float},pseudoCosts::Tuple{Matrix{Float},Matrix{Int}},primalTolerance::Float)::Tuple{Int64,Float}

    bestIndex = 0
    bestScore = 0.0
    for k in 1:length(dscValues)
        score_ = abs(dscValues[k] - round(dscValues[k]))
        score_ = score_*(score_>primalTolerance)

        if score_ > bestScore
            bestIndex = k
            bestScore = score_
        end
    end

    return bestIndex, bestScore
end



function order_of_appearance(dscValues::Vector{Float},pseudoCosts::Tuple{Matrix{Float},Matrix{Int}},primalTolerance::Float;reverse::Bool=false)::Tuple{Int64,Float}

    if reverse
        for k in length(dscValues):-1:1
            if abs(dscValues[k] - round(dscValues[k])) > primalTolerance
                return k, Inf
            end
        end
    else
        for k in 1:length(dscValues)
            if abs(dscValues[k] - round(dscValues[k])) > primalTolerance
                return k, Inf
            end
        end
    end
    return 0, 0.0
end
