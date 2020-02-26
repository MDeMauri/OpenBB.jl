# @Author: Massimo De Mauri <massimo>
# @Date:   2019-05-23T11:47:55+02:00
# @Email:  massimo.demauri@gmail.com
# @Filename: BBnodeChannel.jl
# @Last modified by:   massimo
# @Last modified time: 2020-02-26T17:07:48+01:00
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

    # create a SerialData object pointing to the reserved shared memory space
    serial_ = SerialData(channel.memorySpace)

    # copy the node into the reserved shared memory space
    serialize_in!(serial_,node)

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
    node,_ = BBnode(SerialData(channel.memorySpace))

    # declare the memory space empty
    channel.state[2] = false

    # unlock the memory space
    channel.state[1] = false

    return node
end
