# @Author: Massimo De Mauri <massimo>
# @Date:   2019-03-01T20:22:38+01:00
# @Email:  massimo.demauri@gmail.com
# @Filename: numerical_utilities.jl
# @Last modified by:   massimo
# @Last modified time: 2020-11-06T14:09:43+01:00
# @License: LGPL-3.0
# @Copyright: {{copyright}}

function Infs(n::Integer)::Vector{Float}
    return Inf*ones(n)
end

function Infs(n::Integer,m::Integer)::Matrix{Float}
    return Inf*ones(n,m)
end

function NaNs(n::Integer)::Vector{Float}
    return NaN*ones(n)
end

function NaNs(n::Integer,m::Integer)::Matrix{Float}
    return NaN*ones(n,m)
end


function threshold(number::T1,threshold::T2)::T1 where T1 <: Real where T2 <: Real
    if number >= threshold
        return number
    else
        return T2(0)
    end
end

function clip(number::T1,threshold::T1;relation::Symbol=:<)::T1 where T1 <: Number
    if relation == :< && number < threshold
        return number
    elseif relation == :> && number > threshold
        return relation
    else
        return threshold
    end
end

function clip(number::T1,lower::T1,upper::T1)::T1 where T1 <: Number
    if number < lower
        return lower
    elseif number > upper
        return upper
    else
        return number
    end
end



# cleans an array from the almost zero elements
function cleanZeros!(array::Vector{Float},tolerance::Float)::Nothing
    @. array = array*(abs(array)>tolerance)
    return
end


# perform multiplication assuming Inf*0.0 = 0.0
function weakInfMult(left::Float,right::Float,tolerance::Float=0.0)::Float
    if abs(left) <= tolerance || abs(right) <= tolerance
        return 0.0
    else
        return left*right
    end
end
