# @Author: Massimo De Mauri <massimodemauri>
# @Date:   2019-02-06T16:19:53+01:00
# @Email:  massimo.demauri@gmail.com
# @Project: OpenBB
# @Filename: runtests.jl
# @Last modified by:   massimo
# @Last modified time: 2021-02-12T14:39:16+01:00
# @License: LGPL-3.0
# @Copyright: {{copyright}}

using OpenBB

println("OpenBB tests:")
include("./test_problem_definition_fundamentals.jl");   println(" - problem definition fundamentals, ok")
include("./test_problem_definition_update.jl");         println(" - problem definition update, ok")
include("./test_preprocessing.jl");                     println(" - preprocessing, ok")
include("./test_flat_interface.jl");                    println(" - flat interface, ok")
include("./test_serialization.jl");                     println(" - serialization, ok")
include("./test_LP_subsolvers.jl");                     println(" - LP subsolvers, ok")
include("./test_QP_subsolvers.jl");                     println(" - QP subsolvers, ok")
include("./test_NLP_subsolvers.jl");                    println(" - NLP subsolvers, ok")

if OpenBB.WITH_MPC_ADDON
    include(OpenBB.homeDirectory*"/addons/MPCaddon/test/runtests.jl");  println(" - MPC addon, ok")
end
