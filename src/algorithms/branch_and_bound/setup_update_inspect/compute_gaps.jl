# @Author: Massimo De Mauri <massimo>
# @Date:   2020-01-16T16:18:03+01:00
# @Email:  massimo.demauri@gmail.com
# @Filename: gap_functions.jl
# @Last modified by:   massimo
# @Last modified time: 2021-01-13T12:53:59+01:00
# @License: LGPL-3.0
# @Copyright: {{copyright}}


function compute_gaps(objLoB::Float,objUpB::Float,mode::Symbol=:onUpperBound)::Tuple{Float,Float}
    @assert mode in [:onUpperBound,:onAverage,:onLowerBound]
    if objUpB == Inf || objLoB == -Inf
        return (Inf,Inf)
    elseif mode == :onUpperBound
        absoluteGap = objUpB - objLoB
        relativeGap = absoluteGap/(1e-10 + abs(objUpB))
        return (absoluteGap,relativeGap)
    elseif mode == :onAverage
        absoluteGap = objUpB - objLoB
        relativeGap = absoluteGap/(1e-10 + .5*(abs(objUpB) + abs(objLoB)))
        return (absoluteGap,relativeGap)
    elseif mode == :onLowerBound
        absoluteGap = objUpB - objLoB
        relativeGap = absoluteGap/(1e-10 + abs(objLoB))
        return (absoluteGap,relativeGap)
    end
end

function compute_suboptimality_threshold(objUpB::Float,absGapTol::Float,relGapTol::Float,mode::Symbol=:onUpperBound)::Float
    @assert mode in [:onUpperBound,:onAverage,:onLowerBound]
    if objUpB == Inf
        return Inf
    elseif mode == :onUpperBound
        if objUpB >= 0
            return min(objUpB-absGapTol,(1.0-relGapTol)*objUpB)
        else
            return min(objUpB-absGapTol,(1.0+relGapTol)*objUpB)
        end
    elseif mode == :onAverage
        if objUpB >= 0
            return min(objUpB-absGapTol,(2.0-relGapTol)/(2.0+relGapTol)*objUpB)
        else
            return min(objUpB-absGapTol,(2.0+relGapTol)/(2.0-relGapTol)*objUpB)
        end
    elseif mode == :onLowerBound
        if objLoB >= 0
            return min(objUpB-absGapTol,objUpB/(1.0+relGapTol))
        else
            return min(objUpB-absGapTol,objUpB/(1.0-relGapTol))
        end
    end
end
