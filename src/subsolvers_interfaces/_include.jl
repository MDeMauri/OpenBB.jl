# @Author: Massimo De Mauri <massimo>
# @Date:   2020-01-08T14:39:40+01:00
# @Email:  massimo.demauri@gmail.com
# @Filename: subsolvers_interfaces.jl
# @Last modified by:   massimo
# @Last modified time: 2021-02-12T19:33:32+01:00
# @License: LGPL-3.0
# @Copyright: {{copyright}}

# abstract types
abstract type SubsolverSettings <: AbstractSettings end; struct NullSubsolverSettings <:SubsolverSettings end
abstract type SubsolverWorkspace <: AbstractWorkspace end; struct NullSubsolverWorkspace <: SubsolverWorkspace end

# function that returns the available subsolvers
function get_available_subsolvers()::Array{Vector{String},1}
    out = [["CLP","LP"],["OSQP","LP","QP"],["IPOPT","LP","QP","NLP"]]
    if "QPALM" in keys(installed()) push!(out,["QPALM","LP","QP"]) end
    # if "CPLEX" in keys(installed()) push!(out,["CPLEX","LP","QP"]) end
    # if "Gurobi" in keys(installed()) push!(out,["GUROBI","LP","QP"]) end
    if "JuMP" in keys(installed())
        if "Clp" in keys(installed()) push!(out,["JuMP_CLP","LP"]) end
        if "CPLEX" in keys(installed()) push!(out,["JuMP_CPLEX","LP"]) end
    end
    return out
end

function load_subsolver_interface(name::Symbol)::Nothing
    return load_subsolver_interface(String(name))
end

function load_subsolver_interface(name::String)::Nothing
    if name == "CLP"
        include(homeDirectory*"/src/subsolvers_interfaces/CLP_interface/CLP_interface.jl")
    elseif name[1:4] == "JuMP"
        include(homeDirectory*"/src/subsolvers_interfaces/JuMP_"*uppercase(name[6:end])*"_interface.jl")
    else
        include(homeDirectory*"/src/subsolvers_interfaces/"*uppercase(name)*"_interface.jl")
    end
    return
end


# explaination of the status codes
const statusCodesDictionary = Dict{Int,Symbol}(0=>:optimal,1=>:infeasible,2=>:unreliable,3=>:error)
const reverseStatusCodesDictionary = Dict{Symbol,Int}(:optimal=>0,:infeasible=>1,:unreliable=>2,:error=>3)
