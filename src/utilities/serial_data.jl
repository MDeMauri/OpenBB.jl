# @Author: Massimo De Mauri <massimo>
# @Date:   2020-02-26T17:20:37+01:00
# @Email:  massimo.demauri@gmail.com
# @Filename: serial_data.jl
# @Last modified by:   massimo
# @Last modified time: 2020-12-04T15:12:02+01:00
# @License: LGPL-3.0
# @Copyright: {{copyright}}


# serialization
mutable struct SerialData{T<:AbstractArray}
    data::T
end

function T(serialData::SerialData)::T where T<:AbstractArray
    return T(serialData.data)
end

function Base.lastindex(serial::SerialData)
    return lastindex(serial.data)
end


function Base.getindex(serial::SerialData,index::Int)::Any
    return serial.data[index]
end

function Base.getindex(serial::SerialData,range::UnitRange{Int64})::Any
    return serial.data[range]
end

function Base.setindex!(serial::SerialData,elem::Any,index::Int)::Nothing
    serial.data[index] = elem
    return
end

function Base.setindex!(serial::SerialData{T1},array::T2,range::UnitRange{Int64})::Nothing where T1<:AbstractArray where T2<:AbstractArray
    serial.data[range] = array
    return
end

function Base.length(serial::SerialData)::Int
    return length(serial.data)
end
