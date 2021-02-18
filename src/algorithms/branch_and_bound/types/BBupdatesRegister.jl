# @Author: Massimo De Mauri <massimo>
# @Date:   2019-10-16T22:21:10+02:00
# @Email:  massimo.demauri@gmail.com
# @Filename: BBupdatesRegister.jl
# @Last modified by:   massimo
# @Last modified time: 2020-12-04T14:53:02+01:00
# @License: LGPL-3.0
# @Copyright: {{copyright}}

abstract type AbstractRegister end

mutable struct BBupdatesRegister <: AbstractRegister
    masterVersion::Int64
    updates::Vector{Function}
    arguments::Vector{Tuple}
end


# empty constructor
function BBupdatesRegister()::BBupdatesRegister
    return BBupdatesRegister(0,Function[],Tuple[])
end

# insert new update in the register
function Base.push!(register::BBupdatesRegister,update::Function,arguments::Tuple)::Nothing
    register.masterVersion += 1
    push!(register.updates,update)
    push!(register.arguments,deepcopy(arguments))
    return
end

# take the last update from the register
function Base.pop!(register::BBupdatesRegister)::Tuple{Int64,Function,Tuple}
    register.masterVersion -= 1
    return (register.masterVersion+1,pop!(register.updates),pop!(register.arguments))
end


# eliminates all the updates stored
function reset!(register::BBupdatesRegister)
    register.masterVersion = 0
    empty!(register.updates)
    empty!(register.arguments)
    return
end
