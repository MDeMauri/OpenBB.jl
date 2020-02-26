# @Author: Massimo De Mauri <massimo>
# @Date:   2020-02-26T17:20:37+01:00
# @Email:  massimo.demauri@gmail.com
# @Filename: SerialData.jl
# @Last modified by:   massimo
# @Last modified time: 2020-02-26T18:21:01+01:00
# @License: LGPL-3.0
# @Copyright: {{copyright}}


# serialization
mutable struct SerialData{T<:AbstractArray}
    data::T
end

function T(serialData::SerialData)::T where T<:AbstractArray
    return T(serialData.data)
end

import Base.lastindex
function lastindex(serial::SerialData)
    return lastindex(serial.data)
end


import Base.getindex
function getindex(serial::SerialData,index::Int)::Any
    return serial.data[index]
end

function getindex(serial::SerialData,range::UnitRange{Int64})::Any
    return serial.data[range]
end

import Base.setindex!
function setindex!(serial::SerialData,elem::Any,index::Int)::Nothing
    serial.data[index] = elem
    return
end

function setindex!(serial::SerialData{T1},array::T2,range::UnitRange{Int64})::Nothing where T1<:AbstractArray where T2<:AbstractArray
    serial.data[range] = array
    return
end


import Base.length
function length(serial::SerialData)::Int
    return length(serial.data)
end
