# @Author: Massimo De Mauri <massimodemauri>
# @Date:   2019-02-06T16:19:51+01:00
# @Email:  massimo.demauri@gmail.com
# @Project: OpenBB
# @Filename: OpenBB.jl
# @Last modified by:   massimo
# @Last modified time: 2019-11-22T13:49:31+01:00
# @License: LGPL-3.0
# @Copyright: {{copyright}}


module OpenBB

    info() = print("Hey you, there are no info yet...")

    # include external packages
    using Distributed
    using SparseArrays
    using LinearAlgebra
    using SharedArrays
    using Pkg: installed

    # select the subsolvers to use
    function withOSQP()::Bool
        return "OSQP" in keys(installed())
    end
    function withGUROBI()::Bool
        return "Gurobi" in keys(installed())
    end
    function withQPALM()::Bool
        return "QPALM" in keys(installed())
    end

    function get_available_subsolvers()::Array{String,1}
        out = String[]
        if withOSQP() push!(out,"OSQP") end
        if withQPALM() push!(out,"QPALM") end
        if withGUROBI() push!(out,"GUROBI") end
        return out
    end

    # use or not the MPC addon (the folder containing the mpc toolbox should be placed beside the one containing OpenBB)
    function withMPCaddon()
        return true
    end


    # language interfaces
    include("./problem_definition/problem_definition.jl")

    # solvers
    include("./branch_and_bound/BB.jl")

    # code for preprocessing
    include("./preprocessing/preprocessing.jl")

    # some utilities
    include("./utilities/utilities.jl")

    # subsolvers interfaces
    if withOSQP()
        include("./subsolvers_interfaces/OSQP_interface.jl")
    end
    if withGUROBI()
        include("./subsolvers_interfaces/GUROBI_interface.jl")
    end
    if withQPALM()
        include("./subsolvers_interfaces/QPALM_interface.jl")
    end

    # include heuristics
    include("./heuristics/heuristics.jl")

    # include the flat interface
    include("./alternative_interfaces/flat_interface/flat_interface.jl")

    # load the mpc addon
    if withMPCaddon()
      include(Base.source_path()*"../../../../MPCforOpenBB/src/mpc_addon.jl")
    end

end # module
