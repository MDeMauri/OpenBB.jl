# @Author: Massimo De Mauri <massimodemauri>
# @Date:   2019-02-06T16:19:39+01:00
# @Email:  massimo.demauri@gmail.com
# @Project: OpenBB
# @Filename: BB.jl
# @Last modified by:   massimo
# @Last modified time: 2020-11-24T11:14:48+01:00
# @License: LGPL-3.0
# @Copyright: {{copyright}}


######## include components ######
include("./types/_include.jl")
include("./priority_rules/_include.jl")
include("./pseudo_costs_initialization/_include.jl")
include("./callbacks/_include.jl")
include("./heuristics/_include.jl")
include("./preprocessing/_include.jl")
include("./setup_update_inspect/_include.jl")
include("./solve/_include.jl")
include("./multiprocessing/_include.jl")
