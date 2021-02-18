# @Author: Wim Van Roy
# @Date:   2019-03-11T12:27:34+01:00
# @Email:  wim.vanroy@kuleuven.be
# @Filename: preprocessing.jl
# @Last modified by:   massimo
# @Last modified time: 2020-11-06T18:12:57+01:00
# @License: LGPL-3.0
# @Copyright: {{copyright}}

struct InfeasibleError <: Exception
   msg::String
end

include("./linear_bounds_propagation.jl")
include("./gcd.jl")
include("./sos.jl")
include("./preprocess_root.jl")
include("./preprocess_node.jl")
