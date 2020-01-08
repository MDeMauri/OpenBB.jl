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

    # include the subsolvers
    include("./subsolvers_interfaces/subsolvers_interfaces.jl")

    # include heuristics
    include("./heuristics/heuristics.jl")

    # include the flat interface
    include("./alternative_interfaces/flat_interface/flat_interface.jl")

    # load the mpc addon
    if withMPCaddon()
      include(Base.source_path()*"../../../../MPCforOpenBB/src/mpc_addon.jl")
    end

end # module
