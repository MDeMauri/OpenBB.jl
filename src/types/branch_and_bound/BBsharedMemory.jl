# @Author: Massimo De Mauri <massimo>
# @Date:   2019-05-23T11:37:40+02:00
# @Email:  massimo.demauri@gmail.com
# @Filename: BBsharedMemory.jl
# @Last modified by:   massimo
# @Last modified time: 2019-06-06T15:51:51+02:00
# @License: LGPL-3.0
# @Copyright: {{copyright}}

# abstract and null types
abstract type AbstractSharedMemory end; struct NullSharedMemory <: AbstractSharedMemory end


# channel to send and receive reset_explored_nodes
struct BBnodeChannel <: AbstractChannel{BBnode}
    state::SharedArray{Bool,1}
    memorySpace::SharedArray{Float64,1}
end


function BBnodeChannel(size::Int)
    return BBnodeChannel(SharedArray{Bool,1}([false,false]),SharedArray{Float64,1}(Array{Float64,1}(undef,size)))
end


mutable struct BBsharedMemory{T} <: AbstractSharedMemory where T <: Union{RemoteChannel{Channel{BBnode}},BBnodeChannel}

    inputChannel::T # channel to send nodes
    outputChannel::T # channel to receive nodes
    objectiveBounds::SharedArray{Float64,1} # local objective lower bounds and global objective upper bound
    stats::SharedArray{Int,1} # counters like: how many solutions found
    arrestable::SharedArray{Bool,1} # processes that are done with their local work
end
