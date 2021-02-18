# @Author: Massimo De Mauri <massimo>
# @Date:   2021-01-15T14:17:52+01:00
# @Email:  massimo.demauri@protonmail.com
# @Filename: wra_functions.jl
# @Last modified by:   massimo
# @Last modified time: 2021-02-04T19:48:49+01:00
# @License: LGPL-3.0
# @Copyright: {{copyright}}

# These are the basic building blocks for non-linear constraints and objectives
# Be careful as the output is passed by reference and not copy (for performace reasons)

const Float = Float64
const VirtualVector{T} = Union{Array{T,1},SubArray{T,1}}
const VirtualMatrix{T} = Union{Array{T,2},SubArray{T,2}}

struct WrappedFunction{T<:Number,D}
    baseFunction::Function
    sizeIn::Int
    sizeOut::Union{Tuple{Int},Tuple{Int,Int}}

    # constructors
    function WrappedFunction{T,1}(func::Function,sizeIn::Int,sizeOut::Int)::WrappedFunction{T,1} where T<:Number

        # check input (minimal requirement)
        @assert precompile(func,(Vector{T},Vector{T}))

        # construct helper memory spaces
        input_ = Vector{T}(undef,sizeIn)
        output_ = Vector{T}(undef,sizeOut)


        # wrap the function to check input/output size and to accept SubArrays
        baseFunction =
            if precompile(func,(SubArray{T,1,Array{T,1},Array{Int,1},false},SubArray{T,1,Array{T,1},Tuple{Array{Int,1}},false})) &&
               precompile(func,(SubArray{T,1,Array{T,1},UnitRange{Int},false},SubArray{T,1,Array{T,1},Tuple{Array{Int,1}},false})) &&
               precompile(func,(SubArray{T,1,Array{T,1},Array{Int,1},false},SubArray{T,1,Array{T,1},Tuple{UnitRange{Int}},false})) &&
               precompile(func,(SubArray{T,1,Array{T,1},UnitRange{Int},false},SubArray{T,1,Array{T,1},Tuple{UnitRange{Int}},false}))

                function F1(input::VirtualVector{T},output::VirtualVector{T})::Nothing
                    func(input,output)
                    return
                end

            elseif precompile(func,(Array{T,1},SubArray{T,1,Array{T,1},Tuple{Array{Int,1}},false})) &&
                   precompile(func,(Array{T,1},SubArray{T,1,Array{T,1},Tuple{UnitRange{Int}},false}))

                function F2(input::VirtualVector{T},output::VirtualVector{T})::Nothing
                    @. input_ = input
                    func(input_,output)
                    return
                end
                function F2(input::Vector{T},output::VirtualVector{T})::Nothing
                    func(input,output)
                    return
                end

            elseif precompile(func,(SubArray{T,1,Array{T,1},Array{Int,1},false},Array{T,1})) &&
                   precompile(func,(SubArray{T,1,Array{T,1},UnitRange{Int},false},Array{T,1}))

                function F3(input::VirtualVector{T},output::VirtualVector{T})::Nothing
                    func(input,output_)
                    @. output = output_
                    return
                end
                function F3(input::VirtualVector{T},output::Vector{T})::Nothing
                    func(input,output)
                    return
                end

            else

                function F4(input::VirtualVector{T},output::VirtualVector{T})::Nothing
                    @. input_ = input
                    func(input_,output_)
                    @. output = output_
                    return
                end
                function F4(input::Vector{T},output::VirtualVector{T})::Nothing
                    func(input,output_)
                    @. output = output_
                    return
                end
                function F4(input::VirtualVector{T},output::Vector{T})::Nothing
                    @. input_ = input
                    func(input_,output)
                    return
                end
                function F4(input::Vector{T},output::Vector{T})::Nothing
                    func(input,output)
                    return
                end
            end

        @assert precompile(baseFunction,(SubArray{T,1,Array{T,1},Array{Int,1},false},SubArray{T,1,Array{T,1},Tuple{Array{Int,1}},false})) &&
                precompile(baseFunction,(SubArray{T,1,Array{T,1},UnitRange{Int},false},SubArray{T,1,Array{T,1},Tuple{Array{Int,1}},false})) &&
                precompile(baseFunction,(SubArray{T,1,Array{T,1},Array{Int,1},false},SubArray{T,1,Array{T,1},Tuple{UnitRange{Int}},false})) &&
                precompile(baseFunction,(SubArray{T,1,Array{T,1},UnitRange{Int},false},SubArray{T,1,Array{T,1},Tuple{UnitRange{Int}},false}))
        return new{T,1}(baseFunction,sizeIn,(sizeOut,))

    end
    function WrappedFunction{T,2}(func::Function,sizeIn::Int,sizeOut1::Int,sizeOut2::Int)::WrappedFunction{T,2} where T<:Number

        # check input (minimal requirement)
        @assert precompile(func,(Vector{T},Matrix{T}))

        # define helper memory spaces
        input_ = Vector{T}(undef,sizeIn)
        output_ = Matrix{T}(undef,sizeOut1,sizeOut2)



        # wrap the function to check input/output size and to accept SubArrays
        baseFunction =
            if precompile(func,(SubArray{T,1,Array{T,1},Array{Int,1},false},SubArray{T,2,Array{T,2},Tuple{Array{Int,1},Array{Int,1}},false})) &&
               precompile(func,(SubArray{T,1,Array{T,1},UnitRange{Int},false},SubArray{T,2,Array{T,2},Tuple{Array{Int,1},Array{Int,1}},false})) &&
               precompile(func,(SubArray{T,1,Array{T,1},Array{Int,1},false},SubArray{T,2,Array{T,2},Tuple{UnitRange{Int},UnitRange{Int}},false})) &&
               precompile(func,(SubArray{T,1,Array{T,1},UnitRange{Int},false},SubArray{T,2,Array{T,2},Tuple{UnitRange{Int},UnitRange{Int}},false}))

                function F1(input::VirtualVector{T},output::VirtualMatrix{T})
                    func(input,output)
                    return
                end

            elseif precompile(func,(Array{T,1},SubArray{T,2,Array{T,2},Tuple{Array{Int,1},Array{Int,1}},false})) &&
                   precompile(func,(Array{T,1},SubArray{T,2,Array{T,2},Tuple{UnitRange{Int},UnitRange{Int}},false}))

                function F2(input::VirtualVector{T},output::VirtualMatrix{T})::Nothing
                    @. input_ = input
                    func(input_,output)
                    return
                end
                function F2(input::Vector{T},output::VirtualMatrix{T})::Nothing
                    func(input,output)
                    return
                end

            elseif precompile(func,(SubArray{T,1,Array{T,1},Array{Int,1},false},Array{T,2})) &&
                   precompile(func,(SubArray{T,1,Array{T,1},UnitRange{Int},false},Array{T,2}))


                function F3(input::VirtualVector{T},output::VirtualMatrix{T})::Nothing
                    func(input,output_)
                    @. output = output_
                    return
                end
                function F3(input::VirtualVector{T},output::Matrix{T})::Nothing
                    func(input,output)
                    return
                end

            else
                function F4(input::VirtualVector{T},output::VirtualMatrix{T})::Nothing
                    @. input_ = input
                    func(input_,output_)
                    @. output = output_
                    return
                end
                function F4(input::Vector{T},output::VirtualMatrix{T})::Nothing
                    func(input,output_)
                    @. output = output_
                    return
                end
                function F4(input::VirtualVector{T},output::Matrix{T})::Nothing
                    @. input_ = input
                    func(input_,output)
                    return
                end
                function F4(input::Vector{T},output::Matrix{T})::Nothing
                    func(input,output)
                    return
                end

            end


        @assert precompile(baseFunction,(SubArray{T,1,Array{T,1},Array{Int,1},false},SubArray{T,2,Array{T,2},Tuple{Array{Int,1},Array{Int,1}},false})) &&
                precompile(baseFunction,(SubArray{T,1,Array{T,1},UnitRange{Int},false},SubArray{T,2,Array{T,2},Tuple{Array{Int,1},Array{Int,1}},false})) &&
                precompile(baseFunction,(SubArray{T,1,Array{T,1},Array{Int,1},false},SubArray{T,2,Array{T,2},Tuple{UnitRange{Int},UnitRange{Int}},false})) &&
                precompile(baseFunction,(SubArray{T,1,Array{T,1},UnitRange{Int},false},SubArray{T,2,Array{T,2},Tuple{UnitRange{Int},UnitRange{Int}},false}))
        return new{T,2}(baseFunction,sizeIn,(sizeOut1,sizeOut2))
    end
end

# make the type callable
function (wrapper::WrappedFunction{T,1})(input::VirtualVector{T},output::VirtualVector{T})::Nothing where T<:Number
    @assert length(input) == wrapper.sizeIn
    @assert length(output) == wrapper.sizeOut
    wrapper.baseFunction(input,output)
    return
end
function (wrapper::WrappedFunction{T,2})(input::VirtualVector{T},output::VirtualMatrix{T})::Nothing where T<:Number
    @assert length(input) == wrapper.sizeIn
    @assert size(output) == wrapper.sizeOut
    wrapper.baseFunction(input,output)
    return
end


# copy functions
function Base.deepcopy(wrapper::WrappedFunction{T,D})::WrappedFunction{T,D} where {T<:Number,D}
    return WrappedFunction{T,D}(deepcopy(wrapper.baseFunction),wrapper.sizeIn,wrapper.sizeOut)
end
