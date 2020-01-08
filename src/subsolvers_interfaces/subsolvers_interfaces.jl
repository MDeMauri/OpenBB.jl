
# select the subsolvers interfaces to load
using Pkg: installed
if "Clp" in keys(installed()) include("./CLP_interface.jl") end
if "OSQP" in keys(installed()) include("./OSQP_interface.jl") end
if "QPALM" in keys(installed()) include("./QPALM_interface.jl") end
if "Gurobi" in keys(installed()) include("./GUROBI_interface.jl") end


# function that returns the available subsolvers
function get_available_subsolvers()::Array{String,1}
    out = String[]
    if "Clp" in keys(installed()) push!(out,"CLP") end
    if "OSQP" in keys(installed()) push!(out,"OSQP") end
    if "QPALM" in keys(installed()) push!(out,"QPALM") end
    if "Gurobi" in keys(installed()) push!(out,"GUROBI") end
    return out
end
