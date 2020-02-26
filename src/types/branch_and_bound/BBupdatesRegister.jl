# @Author: Massimo De Mauri <massimo>
# @Date:   2019-10-16T22:21:10+02:00
# @Email:  massimo.demauri@gmail.com
# @Filename: BBupdatesRegister.jl
# @Last modified by:   massimo
# @Last modified time: 2019-10-23T13:16:48+02:00
# @License: LGPL-3.0
# @Copyright: {{copyright}}

abstract type AbstractRegister end

mutable struct BBupdatesRegister <: AbstractRegister
    masterVersion::Int64
    updates::Array{Function,1}
    arguments::Array{Tuple,1}
end


# empty constructor
function BBupdatesRegister()::BBupdatesRegister
    return BBupdatesRegister(0,Function[],Tuple[])
end

# insert new update in the register
import Base.push!
function push!(register::BBupdatesRegister,update::Function,arguments::Tuple)::Nothing
    register.masterVersion += 1
    push!(register.updates,update)
    push!(register.arguments,arguments)
    return
end

# take the last update from the register
import Base.pop!
function pop!(register::BBupdatesRegister)::Tuple{Int64,Function,Tuple}
    register.masterVersion -= 1
    return (register.masterVersion+1,pop!(register.updates),pop!(register.arguments))
end


# eliminates all the updates stored
import Base.empty!
function empty!(register::BBupdatesRegister)
    register.masterVersion = 0
    empty!(register.updates)
    empty!(register.arguments)
    return
end
