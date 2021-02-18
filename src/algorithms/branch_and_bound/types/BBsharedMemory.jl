# @Author: Massimo De Mauri <massimo>
# @Date:   2019-05-23T11:37:40+02:00
# @Email:  massimo.demauri@gmail.com
# @Filename: BBsharedMemory.jl
# @Last modified by:   massimo
# @Last modified time: 2021-01-13T13:29:25+01:00
# @License: LGPL-3.0
# @Copyright: {{copyright}}

# abstract and null types
abstract type AbstractSharedMemory end; struct NullSharedMemory <: AbstractSharedMemory end


# channel to send and receive reset_explored_nodes
struct BBnodeChannel <: AbstractChannel{BBnode}
    state::SharedVector{Bool}
    memorySpace::SharedVector{Float}
end


function BBnodeChannel(size::Int)
    return BBnodeChannel(SharedVector{Bool}([false,false]),SharedVector{Float}(Vector{Float}(undef,size)))
end


mutable struct BBsharedMemory <: AbstractSharedMemory

    inputChannel::BBnodeChannel # channel to send nodes
    outputChannel::BBnodeChannel # channel to receive nodes
    localObjLoBs::SharedVector{Float} # local objective lower bounds
    globalObjUpB::SharedArray{Tuple{Float,Int8},1}   # global objective upper bound and its location
    stats::SharedVector{Int} # counters like: how many solutions found
    arrestable::SharedVector{Bool} # processes that are done with their local work
end


function get_size(channel::BBnodeChannel)::Int
    return length(channel.memorySpace)
end
