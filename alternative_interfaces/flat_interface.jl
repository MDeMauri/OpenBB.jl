# @Author: Massimo De Mauri <massimo>
# @Date:     2019-06-19T16:24:24+02:00
# @Email:    massimo.demauri@gmail.com
# @Filename: flat_interface.jl
# @Last modified by:   massimo
# @Last modified time: 2021-02-09T22:36:01+01:00
# @License: LGPL-3.0
# @Copyright: {{copyright}}

########################################################################################
# this is an addition to the OpenBB interface meant to avoid the use of any struct     #
# in order to facilitate the task of building interfaces for OpenBB in other languages #
########################################################################################
using OpenBB

# define the BBworkspace in the global scope
const FLAT_INTERFACE = true


######################## problem definition ########################
function Problem(problemDict::T)::Problem where T <:Dict
    if isempty(problemDict)
        return NullProblem()
    else
        return Problem(objFun=ObjectiveFunction(problemDict["objFun"]),
                       cnsSet=ConstraintSet(problemDict["cnsSet"]),
                       varSet=VariableSet(problemDict["varSet"]))
    end
end

function ObjectiveFunction(objectiveDict::T)::ObjectiveFunction where T <:Dict
    if isempty(objectiveDict)
        return NullObjective()
    elseif objectiveDict["type"] == "Convex"
        return ConvexObjective(objectiveDict)
    elseif objectiveDict["type"] == "Quadratic"
        return QuadraticObjective(objectiveDict)
    elseif objectiveDict["type"] == "Linear"
        return LinearObjective(objectiveDict)
    else
        error("Objective function type not understood")
    end
end

# dictionary constructors
function LinearObjective(objectiveDict::Dict)::LinearObjective
    @assert objectiveDict["type"] == "Linear"
    if objectiveDict["L"] isa AbstractVector
        L = copy(objectiveDict["L"])
    else
        L = sparsevec(objectiveDict["L"]["inds"],
                      objectiveDict["L"]["vals"],
                      objectiveDict["L"]["n"])
    end
    if "c" in keys(objectiveDict)
        c = objectiveDict["c"]
    else
        c = 0.0
    end
    return OpenBB.LinearObjective(L=L,c=c)
end
function QuadraticObjective(objectiveDict::Dict)::QuadraticObjective
    @assert objectiveDict["type"] == "Quadratic"
    if objectiveDict["Q"] isa AbstractMatrix
        Q = copy(objectiveDict["Q"])
    else
        Q = sparse(objectiveDict["Q"]["rows"],
                   objectiveDict["Q"]["cols"],
                   objectiveDict["Q"]["vals"],
                   objectiveDict["Q"]["n"],
                   objectiveDict["Q"]["m"])
    end
    if objectiveDict["L"] isa AbstractVector
        L = copy(objectiveDict["L"])
    else
        L = sparsevec(objectiveDict["L"]["inds"],
                      objectiveDict["L"]["vals"],
                      objectiveDict["L"]["n"])
    end
    if "c" in keys(objectiveDict)
        c = objectiveDict["c"]
    else
        c = 0.0
    end
    return OpenBB.QuadraticObjective(Q=Q,L=L,c=c)
end
function ConvexObjective(objectiveDict::Dict)::ConvexObjective
    @assert objectiveDict["type"] == "Convex"
    if objectiveDict["grdSparsity"] isa AbstractMatrix
        grdSparsity = copy(objectiveDict["grdSparsity"])
    else
        grdSparsity = sparsevec(objectiveDict["grdSparsity"]["inds"],
                                objectiveDict["grdSparsity"]["vals"],
                                objectiveDict["grdSparsity"]["n"])
    end

    if objectiveDict["hesSparsity"] isa AbstractMatrix
        hesSparsity = copy(objectiveDict["hesSparsity"])
    else
        hesSparsity = sparse(objectiveDict["hesSparsity"]["rows"],
                             objectiveDict["hesSparsity"]["cols"],
                             objectiveDict["hesSparsity"]["vals"],
                             objectiveDict["hesSparsity"]["n"],
                             objectiveDict["hesSparsity"]["m"])
    end

    return OpenBB.ConvexObjective(evalVal=objectiveDict["evalVal"],
                                  evalGrd=objectiveDict["evalGrd"],
                                  evalHes=objectiveDict["evalHes"],
                                  grdSparsity=grdSparsity,hesSparsity=hesSparsity,
                                  typeGrd=objectiveDict["typeGrd"],typeHes=objectiveDict["typeHes"])
end




function ConstraintSet(constraintsDict::T)::ConstraintSet where T <:Dict
    if isempty(constraintsDict)
        return NullConstraintSet()
    elseif constraintsDict["type"] == "Convex"
        return ConvexConstraintSet(constraintsDict)
    elseif constraintsDict["type"] == "Linear"
        return LinearConstraintSet(constraintsDict)
    else
        error("Constraint set type not understood")
    end
end

# dictionary constructors
function LinearConstraintSet(constraintsDict::Dict)::LinearConstraintSet

    @assert constraintsDict["type"] == "Linear"

    if constraintsDict["A"] isa AbstractMatrix
        A = copy(constraintsDict["A"])
    else
        A = sparse(constraintsDict["A"]["rows"],
                   constraintsDict["A"]["cols"],
                   constraintsDict["A"]["vals"],
                   constraintsDict["A"]["n"],
                   constraintsDict["A"]["m"])
    end
    return OpenBB.LinearConstraintSet(A=A,loBs=copy(constraintsDict["loBs"]),upBs=copy(constraintsDict["upBs"]))
end
function ConvexConstraintSet(constraintsDict::Dict)::ConvexConstraintSet

    @assert constraintsDict["type"]=="Convex"

    if constraintsDict["jcbSparsity"] isa AbstractMatrix
        jcbSparsity = copy(constraintsDict["jcbSparsity"])
    else
        jcbSparsity = sparse(constraintsDict["jcbSparsity"]["rows"],
                             constraintsDict["jcbSparsity"]["cols"],
                             constraintsDict["jcbSparsity"]["vals"],
                             constraintsDict["jcbSparsity"]["n"],
                             constraintsDict["jcbSparsity"]["m"])
    end

    numCnss = size(jcbSparsity,1)
    hesSparsity = Vector{SpMatrix{Bool}}(undef,numCnss)
    for k in 1:numCnss
        if constraintsDict["hesSparsity"][k] isa AbstractMatrix
            hesSparsity[k] = copy(constraintsDict["hesSparsity"][k])
        else
            hesSparsity[k] = sparse(constraintsDict["hesSparsity"][k]["rows"],
                                    constraintsDict["hesSparsity"][k]["cols"],
                                    constraintsDict["hesSparsity"][k]["vals"],
                                    constraintsDict["hesSparsity"][k]["n"],
                                    constraintsDict["hesSparsity"][k]["m"])
        end
    end

    return OpenBB.ConvexConstraintSet(evalVal=constraintsDict["evalVal"],
                                      evalJcb=constraintsDict["evalJcb"],
                                      evalHes=constraintsDict["evalHes"],
                                      jcbSparsity=jcbSparsity,
                                      hesSparsity=hesSparsity,
                                      typeJcb=constraintsDict["typeJcb"],
                                      typeHes=constraintsDict["typeHes"],
                                      loBs=constraintsDict["loBs"],
                                      upBs=constraintsDict["upBs"])
end





function VariableSet(variableDict::T)::VariableSet where T <:Dict

    if isempty(variableDict)
        return EmptyVariableSet()
    else
        if "vals" in keys(variableDict)
            vals = copy(variableDict["vals"])
        else
            vals = Float64[]
        end
        if "dscIndices" in keys(variableDict)
            dscIndices = variableDict["dscIndices"]
        else
            dscIndices = Int[]
        end
        if "sos1Groups" in keys(variableDict)
            sos1Groups = copy(variableDict["sos1Groups"])
        else
            sos1Groups = Int.(zeros(length(dscIndices)))
        end
        if "pseudoCosts" in keys(variableDict)
            pseudoCosts = copy(variableDict["pseudoCosts"])
        else
            pseudoCosts = (1e4.*ones(length(dscIndices),2),Int.(zeros(length(dscIndices),2)))
        end
        return VariableSet(loBs=copy(variableDict["loBs"]),upBs=copy(variableDict["upBs"]),vals=vals,
                                     dscIndices=dscIndices,sos1Groups=sos1Groups,pseudoCosts=pseudoCosts)
    end
end


######################## setup ########################
# shorthand function to collect subsolver settings
function SubsolverSettings(settingsDict::T)::SubsolverSettings where T <:Dict
    # convert settings for subsolver
    load_subsolver_interface(settingsDict["subsolverName"])
    ssSettings = OpenBB.eval(Symbol(settingsDict["subsolverName"],"settings"))()
    for field in fieldnames(typeof(ssSettings))
        if string(field) in keys(settingsDict)
            setfield!(ssSettings,field,settingsDict[string(field)])
        end
    end
    return ssSettings
end

# shorthand function to collect BB heuristic settings
function BBheuristicsSettings(settingsDict::T)::BBheuristicsSettings where T <:Dict
    bbHSettings = eval(Symbol("BB",settingsDict["heuristicsName"],"Settings"))
    for field in fieldnames(typeof(bbHSettings))
        if string(field) in keys(settingsDict)
            setfield!(bbHSettings,field,settingsDict[string(field)])
        end
    end
    return bbHSettings
end


# shorthand function to collect BB settings
function BBsettings(settingsDict::T)::BBsettings where T <:Dict

    # default BB settings
    bbSettings = BBsettings()

    # collect user settings
    for field in fieldnames(BBsettings)
        if string(field) in keys(settingsDict)
            if field == :subsolverSettings
                bbSettings.subsolverSettings = SubsolverSettings(settingsDict["subsolverSettings"])
            elseif field == :heuristicsSettings
                bbSettings.heuristicsSettings = BBheuristicsSettings(string(field))
            elseif fieldtype(BBsettings,field) <: Symbol
                setfield!(bbSettings,field,Symbol(settingsDict[string(field)]))
            elseif fieldtype(BBsettings,field) <: Tuple{Symbol,Vararg}
                setfield!(bbSettings,field,tuple(Symbol(settingsDict[string(field)][1],settingsDict[string(field)][2:end]...)))
            else
                setfield!(bbSettings,field,settingsDict[string(field)])
            end
        end
    end

    return bbSettings
end

# shorthand function to collect HBB heuristic settings
function HBBheuristicsSettings(settingsDict::T)::HBBheuristicsSettings where T <:Dict
    hbbHSettings = eval(Symbol("HBB",settingsDict["heuristicsName"],"Settings"))
    for field in fieldnames(typeof(hbbHSettings))
        if string(field) in keys(settingsDict)
            setfield!(hbbHSettings,field,settingsDict[string(field)])
        end
    end
    return hbbHSettings
end



# shorthand function to collect HBB settings
function HBBsettings(settingsDict::T)::HBBsettings where T <:Dict

    # default HBB settings
    hbbSettings = HBBsettings()

    # collect user settings
    for field in fieldnames(HBBsettings)
        if string(field) in keys(settingsDict)
            if field == :mipSettings
                hbbSettings.mipSettings=BBsettings(settingsDict["mipSettings"])
            elseif field == :nlpSettings
                hbbSettings.nlpSettings=SubsolverSettings(settingsDict["nlpSettings"])
            elseif field == :heuristicSettings
                hbbSettings.heuristicSettings=HBBheuristicsSettings(settingsDict["nlpSettings"])
            elseif fieldtype(HBBsettings,field) <: Symbol
                setfield!(hbbSettings,field,Symbol(settingsDict[string(field)]))
            elseif fieldtype(HBBsettings,field) <: Tuple{Symbol,Vararg}
                setfield!(hbbSettings,field,tuple(Symbol(settingsDict[string(field)][1],settingsDict[string(field)][2:end]...)))
            else
                setfield!(hbbSettings,field,settingsDict[string(field)])
            end
        end
    end

    return hbbSettings
end



# setup a workspace
function setup(algorithm::String,problemDict::T1,settingsDict::T2)::SupersolverWorkspace where T1 <:Dict where T2 <:Dict

    if algorithm == "BB" || algorithm == "bb"
        workspace = setup(Problem(problemDict),BBsettings(settingsDict))
    elseif algorithm == "HBB" || algorithm == "hbb"
        workspace = setup(Problem(problemDict),HBBsettings(settingsDict))
    else
        error("setup: algorithm unknown")
    end

    return workspace
end


######################## solve ########################
function solveB(workspace::SupersolverWorkspace,initialPoint::Vector{Float}=Float[])::Nothing
    if workspace isa BBworkspace
        solve!(workspace)
    else
        solve!(workspace,initialPoint)
    end
    return
end




######################## inspect ########################

# ...
function get_status_dict(workspace::SupersolverWorkspace)::Dict{String,Any}
    status = get_status(workspace)
    out = Dict{String,Any}()
    for fn in fieldnames(typeof(status))
        out[String(fn)] = getfield(status,fn)
    end
    return out
end

function get_best_feasible_node_dict(workspace::SupersolverWorkspace)::Dict{String,Any}
    node = get_best_feasible_node(workspace)
    out = Dict{String,Any}()
    if !(node isa NullBBnode)
        for field in fieldnames(BBnode)
            out[String(field)] = getfield(node,field)
        end
    end
    return out
end


function get_all_feasible_nodes_dict(workspace::SupersolverWorkspace)::Array{Dict{String,Any},1}
    tmp = get_all_feasible_nodes(workspace)
    out = Array{Dict{String,Any},1}(undef,length(tmp))
    for k in length(tmp)
        out[k] = Dict{String,Any}()
        for field in fieldnames(BBnode)
            out[k][String(field)] = getfield(tmp[k],field)
        end
    end
    return out
end

function get_best_node_dict(workspace::SupersolverWorkspace)::Dict{String,Any}
    if workspace isa HBBworkspace
        return get_best_feasible_node_dict(workspace)
    end
    node = get_best_node(workspace)
    out = Dict{String,Any}()
    if !(node isa NullBBnode)
        for field in fieldnames(BBnode)
            out[String(field)] = getfield(node,field)
        end
    end
    return out
end



function get_constraints_dict(workspace::SupersolverWorkspace)::Dict{String,Any}
    tmp = get_constraints(workspace)
    out = Dict{String,Any}()
    for field in fieldnames(typeof(tmp))
        out[String(field)] = getfield(tmp,field)
    end
    return out
end

function get_objective_dict(workspace::SupersolverWorkspace)::Dict{String,Any}

    tmp = get_objective(workspace)
    out = Dict{String,Any}()
    for field in fieldnames(typeof(tmp))
        out[String(field)] = getfield(tmp,field)
    end
    return out
end

# ...


######################## update workspace ########################

function updateB(workspace::SupersolverWorkspace)::Nothing
    return update!(workspace)
end

function resetB(workspace::SupersolverWorkspace)::Nothing
    return reset!(workspace)
end

function clearB(workspace::SupersolverWorkspace)::Nothing
    return clear!(workspace)
end




######################## update problem ########################

function append_constraintsB(workspace::SupersolverWorkspace,constraintsDict::T,suppressErrors::Bool=false)::Nothing where T <:Dict
    return append_constraints!(workspace,ConstraintSet(constraintsDict),suppressErrors=suppressErrors)
end

function insert_constraintsB(workspace::SupersolverWorkspace,constraintsDict::T,index::Int,suppressErrors::Bool=false)::Nothing where T <:Dict
    return insert_constraints!(workspace,index,ConstraintSet(constraintsDict),suppressErrors=suppressErrors)
end

function remove_constraintsB(workspace::SupersolverWorkspace,indices::Vector{Int},suppressErrors::Bool=false)::Nothing
    return remove_constraints!(workspace,indices,suppressErrors=suppressErrors)
end


function permute_constraintsB(workspace::SupersolverWorkspace,permutation::Vector{Int},suppressErrors::Bool=false)::Nothing
    return permute_constraints!(workspace,permutation,suppressErrors=suppressErrors)
end

function update_boundsB(workspace::SupersolverWorkspace,boundsDict::T,suppressErrors::Bool=false)::Nothing where T <:Dict
    boundsIn = Dict{Symbol,Vector{Float64}}()
    for key in keys(boundsDict)
        boundsIn[Symbol(key)] = boundsDict[key]
    end
    return update_bounds!(workspace; boundsIn...,suppressErrors=suppressErrors)
end


function set_objectiveB(workspace::SupersolverWorkspace,newObjectiveDict::T,suppressErrors::Bool=false)::Nothing where T <:Dict

    return set_objective!(workspace,ObjectiveFunction(newObjectiveDict),suppressErrors=suppressErrors)
end

function set_constraintSetB(workspace::SupersolverWorkspace,newConstraintSetDict::T,suppressErrors::Bool=false)::Nothing where T <:Dict

    return set_constraintSet!(workspace,ConstraintSet(newConstraintSetDict),suppressErrors=suppressErrors)
end


function append_problemB(workspace::SupersolverWorkspace,problemDict::T,suppressErrors::Bool=false)::Nothing where T <:Dict
    return append_problem!(workspace,Problem(problemDict),suppressErrors=suppressErrors)
end


function integralize_variablesB(workspace::SupersolverWorkspace,newDscIndices::Vector{Int},newSos1Groups::Vector{Int}=Int[],suppressErrors::Bool=false)::Nothing
    return integralize_variables!(workspace,newDscIndices,suppressErrors=suppressErrors)
end

######################## update settings ########################
function update_objectiveCutoffB(workspace::SupersolverWorkspace,newCutoff::Float64,suppressErrors::Bool=false)::Nothing
    return update_objectiveCutoff!(workspace,newCutoff,suppressErrors=suppressErrors)
end


function addto_blackListB(workspace::SupersolverWorkspace,assignment::Vector{Float64})
    return addto_blackList!(workspace,assignment)
end



######################## optimal control mode updates ########################
function sort_constraints_ocB(workspace::SupersolverWorkspace,suppressErrors::Bool=false)::Nothing
    return sort_constraints_oc!(workspace,suppressErrors=suppressErrors)
end

function insert_constraints_ocB(workspace::SupersolverWorkspace,newConstraintsDict::Dict,suppressErrors::Bool=false)::Nothing
    return insert_constraints_oc!(workspace,ConstraintSet(newConstraintsDict),suppressErrors=suppressErrors)
end
