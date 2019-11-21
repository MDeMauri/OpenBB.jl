# @Author: Massimo De Mauri <massimo>
# @Date:   2019-05-23T11:47:55+02:00
# @Email:  massimo.demauri@gmail.com
# @Filename: communication.jl
# @Last modified by:   massimo
# @Last modified time: 2019-11-21T14:13:19+01:00
# @License: LGPL-3.0
# @Copyright: {{copyright}}


# check if the required location of the channel is locked by another process
import Base.isready
function isready(channel::BBnodeChannel)::Bool
    if channel.state[2]
        return true
    else
        return false
    end
end


# place a node on the shared memory
import Base.put!
function put!(channel::BBnodeChannel,node::T;timeout::Float64=Inf)::Nothing where T <: AbstractBBnode

    # wait for the channel to be empty and unlocked
    startTime = time()
    while true
        if time() - startTime > timeout
            @error "communication timeout"
        end
        if !channel.state[1] && !isready(channel)
            break
        end
        s = time()
        while time() - s < 0.0001
        end
    end

    # lock the memory space
    channel.state[1] = true

    # copy the node into the reserved shared memory space
    flatten_in!(node,channel.memorySpace)

    # declare the memory space full
    channel.state[2] = true

    # unlock the memory space
    channel.state[1] = false

    return
end


import Base.take!
function take!(channel::BBnodeChannel;timeout::Float64=Inf)

    # wait for the process to be full and unlocked
    startTime = time()
    while true
        if time()-startTime>timeout
            @error "communication timeout"
        end
        if  !channel.state[1] && isready(channel)
            break
        end
        s = time()
        while time() - s < 0.0001
        end
    end

    # lock the memory space
    channel.state[1] = true

    # read the memory space
    node = rebuild_node(channel.memorySpace)

    # declare the memory space empty
    channel.state[2] = false

    # unlock the memory space
    channel.state[1] = false

    return node
end
