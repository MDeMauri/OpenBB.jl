# @Author: Massimo De Mauri <massimo>
# @Date:   2019-06-07T18:58:51+02:00
# @Email:  massimo.demauri@gmail.com
# @Filename: pseudo_costs_initialization.jl
# @Last modified by:   massimo
# @Last modified time: 2021-01-08T19:27:55+01:00
# @License: LGPL-3.0
# @Copyright: {{copyright}}

# wrappers for initialization functions
function initialize_pseudoCosts!(functionTuple::Tuple{Symbol,Vararg},pseudoCosts::Tuple{Matrix{Float},Matrix{Int}})::Nothing
    return OpenBB.eval(functionTuple[1])(pseudoCosts,functionTuple[2:end]...)
end


# actual initialization functions
function initialize_to_constant!(pseudoCosts::Tuple{Matrix{Float},Matrix{Int}},value::Float)::Nothing
    @. pseudoCosts[1] = value
    @. pseudoCosts[2] = 0
    return
end
