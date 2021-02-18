# @Author: Massimo De Mauri <massimo>
# @Date:   2019-05-23T11:47:55+02:00
# @Email:  massimo.demauri@gmail.com
# @Filename: BBnodeChannel.jl
# @Last modified by:   massimo
# @Last modified time: 2020-12-11T16:18:37+01:00
# @License: LGPL-3.0
# @Copyright: {{copyright}}




# check if the required location of the channel is locked by another process
function Base.islocked(channel::BBnodeChannel)::Bool
    return channel.state[1]
end
function Base.isempty(channel::BBnodeChannel)::Bool
    return !channel.state[2]
end
function Base.isready(channel::BBnodeChannel)::Bool
    return !islocked(channel) && !isempty(channel)
end

# place a node on the shared memory
function Base.put!(channel::BBnodeChannel,node::T;timeout::Float=Inf,mode::Symbol=:normal)::Nothing where T <: AbstractBBnode

    # wait for the channel to be empty and unlocked
    startTime = time()
    while islocked(channel) || (!isempty(channel) && mode!=:forced)
        if time() - startTime > timeout
            error("communication timeout")
        else
            sleep(0.001)
        end
    end

    # lock the memory space
    channel.state[1] = true

    # declare the memory space full
    channel.state[2] = true

    # create a SerialData object pointing to the reserved shared memory space
    serial_ = SerialData(channel.memorySpace)

    # copy the node into the reserved shared memory space
    serialize_in!(serial_,node)

    # unlock the memory space
    channel.state[1] = false

    return
end


function Base.take!(channel::BBnodeChannel;timeout::Float=Inf)

    # wait for the process to be full and unlocked
    startTime = time()
    while islocked(channel) || isempty(channel)
        if time()-startTime>timeout
            error("communication timeout")
        else
            sleep(0.001)
        end
    end

    # lock the memory space
    channel.state[1] = true

    # declare the memory space empty
    channel.state[2] = false

    # read the memory space
    node,_ = BBnode(SerialData(channel.memorySpace))

    # unlock the memory space
    channel.state[1] = false

    return node
end
