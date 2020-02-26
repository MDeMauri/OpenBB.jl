# @Author: Massimo De Mauri <massimodemauri>
# @Date:   2019-02-06T16:19:39+01:00
# @Email:  massimo.demauri@gmail.com
# @Project: OpenBB
# @Filename: BB.jl
# @Last modified by:   massimo
# @Last modified time: 2020-02-26T21:29:12+01:00
# @License: LGPL-3.0
# @Copyright: {{copyright}}


######## include components ######
include("./priority_rules/priority_rules.jl")
include("./pseudo_costs_initialization/pseudo_costs_initialization.jl")
include("./setup_update_inspect/setup_update_inspect.jl")
include("./solve/solve!.jl")
include("./multiprocessing/multiprocessing.jl")
