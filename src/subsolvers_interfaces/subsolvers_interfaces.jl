# @Author: Massimo De Mauri <massimo>
# @Date:   2020-01-08T14:39:40+01:00
# @Email:  massimo.demauri@gmail.com
# @Filename: subsolvers_interfaces.jl
# @Last modified by:   massimo
# @Last modified time: 2020-01-08T16:39:56+01:00
# @License: LGPL-3.0
# @Copyright: {{copyright}}


# select the subsolvers interfaces to load
using Pkg: installed
if "Clp" in keys(installed()) include("./CLP_interface.jl") end
if "OSQP" in keys(installed()) include("./OSQP_interface.jl") end
if "QPALM" in keys(installed()) include("./QPALM_interface.jl") end
if "Gurobi" in keys(installed()) include("./GUROBI_interface.jl") end


# function that returns the available subsolvers
function get_available_subsolvers()::Array{Array{String,1},1}
    out = Array{String,1}[]
    if "Clp" in keys(installed()) push!(out,["CLP","LP"]) end
    if "OSQP" in keys(installed()) push!(out,["OSQP","QP"]) end
    if "QPALM" in keys(installed()) push!(out,["QPALM","QP"]) end
    if "Gurobi" in keys(installed()) push!(out,["GUROBI","QP"]) end
    return out
end
