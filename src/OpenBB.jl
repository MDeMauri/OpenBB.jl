# @Author: Massimo De Mauri <massimodemauri>
# @Date:   2019-02-06T16:19:51+01:00
# @Email:  massimo.demauri@gmail.com
# @Project: OpenBB
# @Filename: OpenBB.jl
# @Last modified by:   massimo
# @Last modified time: 2021-02-12T18:14:25+01:00
# @License: LGPL-3.0
# @Copyright: {{copyright}}

__precompile__()

module OpenBB

    info() = print("Hey you, there are no info yet...")

    # include external packages
    using Distributed
    using SparseArrays
    using LinearAlgebra
    using Statistics: median, mean
    using SharedArrays
    using Pkg: dependencies
    using Base: RefValue

    # find src folder
    const homeDirectory = Base.source_dir()[1:end-4]
    const tempDirectory = homeDirectory*"/non_public_temp_bin/"

    # container for error logging
    errorLogs = Vector{String}()
    function get_errorLogs()::Nothing
        @everywhere try @. println(OpenBB.errorLogs) catch err end
        return
    end

    # type aliases
    const Float = Float64
    const SpMatrix{T} = SparseMatrixCSC{T,Int}
    const SpVector{T} = SparseVector{T,Int}
    const VirtualVector{T} = Union{Array{T,1},SubArray{T,1}}
    const VirtualMatrix{T} = Union{Array{T,2},SubArray{T,2}}
    const AbMatrix{T} = Union{Matrix{T},SpMatrix{T}}
    const AbVector{T} = Union{Vector{T},SpVector{T}}


    # define infinity for ints
    const InfI64 = 2^63-1
    const InfI32 = 2^31-1

    # base abstract types
    abstract type AbstractSettings end; struct NullSettings <: AbstractSettings end
    abstract type AbstractWorkspace end; struct NullWorkspace <: AbstractWorkspace end
    abstract type AbstractNode end; struct NullNode <: AbstractNode end

    # some utilities
    include(homeDirectory*"/src/utilities/_include.jl")

    # custom types
    include(homeDirectory*"/src/problem_definition/_include.jl")

    # include the subsolvers
    include(homeDirectory*"/src/subsolvers_interfaces/_include.jl")

    # Branch&Bound and Hybrid-Branch&Bound routines
    include(homeDirectory*"/src/algorithms/_include.jl")

    # include addons
    include(homeDirectory*"/addons/_include.jl")

    # include default subsolvers:
    load_subsolver_interface("CLP")
    load_subsolver_interface("OSQP")
    load_subsolver_interface("IPOPT")

    # optional inclusions
    function include_flat_interface()::Nothing
        include(homeDirectory*"/alternative_interfaces/flat_interface.jl")
        if WITH_MPC_ADDON
            include(homeDirectory*"/addons/MPCaddon/alternative_interfaces/flat_interface.jl")
        end
        return
    end

    # test function
    function runtests()::Nothing
        Main.include(homeDirectory*"/test/runtests.jl")
        return
    end

end # module
