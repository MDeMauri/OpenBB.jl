# @Author: Massimo De Mauri <massimo>
# @Date:   2019-06-19T16:24:24+02:00
# @Email:  massimo.demauri@gmail.com
# @Filename: flat_interface.jl
# @Last modified by:   massimo
# @Last modified time: 2019-06-20T00:55:10+02:00
# @License: LGPL-3.0
# @Copyright: {{copyright}}

############################################################################
# this is an interface which does not use any struct in order to           #
# facilitate the task of building interfaces for OpenBB in other languages #
############################################################################


using OpenBB

# define the BBworkspace in the global scope
global globalWorkspace = NullWorkspace()

######################## problem definition ########################

function ObjectiveFunction(objectiveDict::Dict)::AbstractObjectiveFunction
  if isempty(objectiveDict)
  return NullObjectiveFunction()
  elseif "Q" in keys(objectiveDict) && "L" in keys(objectiveDict)
  return QuadraticObjective(Q=copy(objectiveDict["Q"]),L=copy(objectiveDict["L"]))
  elseif "L" in keys(objectiveDict)
  return LinearObjective(L=copy(objectiveDict["L"]))
  else
    @error "Objective function type not understood"
  return NullObjectiveFunction()
  end
end


function ConstraintSet(constraintsDict::Dict)::AbstractConstraintSet
  if isempty(constraintsDict)
    return NullConstraintSet()
  elseif "A" in keys(constraintsDict) && "loBs" in keys(constraintsDict) && "upBs" in keys(constraintsDict)
      return LinearConstraintSet(A=copy(constraintsDict["A"]),
                                        loBs=copy(constraintsDict["loBs"]),
                                        upBs=copy(constraintsDict["upBs"]))
  else
    @error "Constraint set type not understood"
    return NullConstraintSet()
  end
end

function VariableSet(variableDict::Dict)::AbstractVariableSet

  if isempty(variableDict)
    return NullVariableSet()
  else

    if "dscIndices" in keys(variableDict)
      dscIndices = copy(variableDict["dscIndices"])
    else
      dscIndices = Int[]
    end
    if "sos1Groups" in keys(variableDict)
      sos1Groups = copy(variableDict["sos1Groups"])
    else
      sos1Groups = Int[]
    end
    if "pseudoCosts" in keys(variableDict)
      pseudoCosts = copy(variableDict["pseudoCosts"])
    else
      pseudoCosts = Array{Float64,2}(undef,0,2)
    end
    return VariableSet(loBs=copy(variableDict["loBs"]),upBs=copy(variableDict["upBs"]),vals=copy(variableDict["vals"]),
                              dscIndices=dscIndices,sos1Groups=sos1Groups,pseudoCosts=pseudoCosts)
  end
end



function Problem(problemDict::Dict)::Problem
  if isempty(problemDict)
    return NullProblem()
  else
    return Problem(objFun=ObjectiveFunction(problemDict["objFun"]),
                          cnsSet=ConstraintSet(problemDict["cnsSet"]),
                          varSet=VariableSet(problemDict["varSet"]))
  end
end


######################## setup ########################

# setup a globalWorkspace
function setup(subsolver::String,problemDict::Dict,bbSettingsDict::Dict,ssSettingsDict::Dict)::Nothing

  # manage settings for BB algorithm
  bb_settings = BBsettings()
  for pair in bbSettingsDict
    setfield!(bb_settings,Symbol(pair[1]),pair[2])
  end

  # manage subsolver settings
  if subsolver == "osqp"
    ss_settings = OSQPsettings()
  elseif subsolver == "gurobi"
    ss_settings = GUROBIsettings()
  elseif subsolver == "qpalm"
    ss_settings = QPALMsettings()
  else
    @error "Unknown subsolver"
  end

  for pair in ssSettingsDict
    setfield!(ss_settings,Symbol(pair[1]),pair[2])
  end

  # create a OpenBB problem description
  global globalWorkspace = setup(Problem(problemDict),bb_settings,ss_settings)

  return
end





######################## solve ########################

function solve_b()::Nothing
  if globalWorkspace isa NullWorkspace
    @error "workspace not initialized, please run setup(problemDict,bbSettingsDict,ssSettingsDict)"
  end
  solve!(globalWorkspace)
  return
end





######################## inspect ########################

# ...
function get_settings()::Dict{String,Any}
  tmp = get_settings(globalWorkspace)
  out = Dict{String,Any}()
  for field in fieldnames(BBsettings)
    out[String(field)] = getfield(tmp,field)
  end
  return out
end

# ...
function get_subsolver_settings()::Dict{String,Any}
  tmp = get_settings(globalWorkspace)
  out = Dict{String,Any}()
  for field in fieldnames(typeof(tmp))
    out[String(field)] = getfield(tmp,field)
  end
  return out
end

# ...
function get_status()::Dict{String,Any}
  if globalWorkspace isa NullWorkspace
    @error "workspace not initialized, please run setup(problemDict,bbSettingsDict,ssSettingsDict)"
  end
  status = get_status(globalWorkspace)
  out = Dict{String,Any}()
  for fn in fieldnames(BBstatus)
    out[String(fn)] = getfield(status,fn)
  end
  return out
end

function get_best_solution(localOnly::Bool=false)::Dict{String,Any}
  if globalWorkspace isa NullWorkspace
    @error "workspace not initialized, please run setup(problemDict,bbSettingsDict,ssSettingsDict)"
  end
  node = get_best_solution(globalWorkspace,localOnly=localOnly)
  out = Dict{String,Any}()
  if !(node isa NullBBnode)
    for field in fieldnames(BBnode)
      out[String(field)] = getfield(node,field)
    end
  end
  return out
end


function get_all_solutions(localOnly::Bool=false)::Array{Dict{String,Any},1}
  if globalWorkspace isa NullWorkspace
    @error "workspace not initialized, please run setup(problemDict,bbSettingsDict,ssSettingsDict)"
  end
  tmp = get_all_solutions(globalWorkspace,localOnly=localOnly)
  out = Array{Dict{String,Any},1}(undef,length(tmp))
  for k in length(tmp)
      out[k] = Dict{String,Any}()
      for field in fieldnames(BBnode)
        out[k][String(field)] = getfield(tmp[k],field)
      end
  end
  return out
end

function get_best_node(localOnly::Bool=false)::Dict{String,Any}
  if globalWorkspace isa NullWorkspace
    @error "workspace not initialized, please run setup(problemDict,bbSettingsDict,ssSettingsDict)"
  end
  node = get_best_node(globalWorkspace,localOnly=localOnly)
  out = Dict{String,Any}()
  if !(node isa NullBBnode)
    for field in fieldnames(BBnode)
      out[String(field)] = getfield(node,field)
    end
  end
  return out


end

function get_numVariables()::Int
  if globalWorkspace isa NullWorkspace
    @error "workspace not initialized, please run setup(problemDict,bbSettingsDict,ssSettingsDict)"
  end
  return get_numVariables(globalWorkspace)
end

function get_numDiscreteVariables()::Int
  if globalWorkspace isa NullWorkspace
    @error "workspace not initialized, please run setup(problemDict,bbSettingsDict,ssSettingsDict)"
  end
  return get_numDiscreteVariables(globalWorkspace)
end

function get_numConstraints()::Int
  if globalWorkspace isa NullWorkspace
    @error "workspace not initialized, please run setup(problemDict,bbSettingsDict,ssSettingsDict)"
  end
  return get_numConstraints(globalWorkspace)
end

function get_constraints_sparsity()::Tuple{Array{Int,1},Array{Int,1}}
  if globalWorkspace isa NullWorkspace
    @error "workspace not initialized, please run setup(problemDict,bbSettingsDict,ssSettingsDict)"
  end
  return get_constraints_sparsity(globalWorkspace)
end

function get_constraint_sparsity(index::Int)::Array{Int,1}
  if globalWorkspace isa NullWorkspace
    @error "workspace not initialized, please run setup(problemDict,bbSettingsDict,ssSettingsDict)"
  end
  return get_constraint_sparsity(globalWorkspace,index)
end

function get_objective_sparsity()::Tuple{Array{Int,1},Array{Int,1}}
  if globalWorkspace isa NullWorkspace
    @error "workspace not initialized, please run setup(problemDict,bbSettingsDict,ssSettingsDict)"
  end
  return get_objective_sparsity(globalWorkspace)
end

function get_variableBounds()::Tuple{Array{Float64,1},Array{Float64,1}}
  if globalWorkspace isa NullWorkspace
    @error "workspace not initialized, please run setup(problemDict,bbSettingsDict,ssSettingsDict)"
  end
  return get_variableBounds(globalWorkspace)
end

function get_constraintBounds()::Tuple{Array{Float64,1},Array{Float64,1}}
  if globalWorkspace isa NullWorkspace
    @error "workspace not initialized, please run setup(problemDict,bbSettingsDict,ssSettingsDict)"
  end
  return get_constraintBounds(globalWorkspace)
end

# ...
function get_numActiveNodes()::Int
  return length(globalWorkspace.activeQueue)
end

# ...
function get_numUnactiveNodes()::Int
  return length(globalWorkspace.unactivePool)
end

# ...
function get_numSolutions()::Int
  return length(globalWorkspace.solutionPool)
end





######################## update workspace ########################

function reset_explored_nodes_b(localOnly::Bool=false)::Nothing
  if globalWorkspace isa NullWorkspace
    @error "workspace not initialized, please run setup(problemDict,bbSettingsDict,ssSettingsDict)"
  end
  return reset_explored_nodes!(globalWorkspace,localOnly=localOnly)
end

function update_b(localOnly::Bool=false)::Nothing
  if globalWorkspace isa NullWorkspace
    @error "workspace not initialized, please run setup(problemDict,bbSettingsDict,ssSettingsDict)"
  end
  return update!(globalWorkspace,localOnly=localOnly)
end

function reset_b(localOnly::Bool=false)::Nothing
  if globalWorkspace isa NullWorkspace
    @error "workspace not initialized, please run setup(problemDict,bbSettingsDict,ssSettingsDict)"
  end
  return reset!(globalWorkspace,localOnly=localOnly)
end

function clear_b(localOnly::Bool=false)::Nothing
  if globalWorkspace isa NullWorkspace
    @error "workspace not initialized, please run setup(problemDict,bbSettingsDict,ssSettingsDict)"
  end
  return clear!(globalWorkspace,localOnly=localOnly)
end




######################## update problem ########################

function append_constraints_b(constraintsDict::Dict,
                              suppressWarnings::Bool=false,
                              suppressUpdate::Bool=false,
                              localOnly::Bool=false)::Nothing
  if globalWorkspace isa NullWorkspace
    @error "workspace not initialized, please run setup(problemDict,bbSettingsDict,ssSettingsDict)"
  end
  return append_constraints!(globalWorkspace,constraintsDict["A"],constraintsDict["loBs"],constraintsDict["upBs"],
                                    suppressWarnings=suppressWarnings,
                                    suppressUpdate=suppressUpdate,
                                    localOnly=localOnly)
end

function insert_constraints_b(constraintsDict::Dict,index::Int,
                              suppressWarnings::Bool=false,
                              suppressUpdate::Bool=false,
                              localOnly::Bool=false)::Nothing
  if globalWorkspace isa NullWorkspace
    @error "workspace not initialized, please run setup(problemDict,bbSettingsDict,ssSettingsDict)"
  end
  return insert_constraints!(globalWorkspace,index,constraintsDict["A"],constraintsDict["cnsLoBs"],constraintsDict["cnsUpBs"],
                                    suppressWarnings=suppressWarnings,
                                    suppressUpdate=suppressUpdate,
                                    localOnly=localOnly)
end

function remove_constraints_b(indices::Array{Int,1},
                              suppressWarnings::Bool=false,
                              suppressUpdate::Bool=false,
                              localOnly::Bool=false)::Nothing
  if globalWorkspace isa NullWorkspace
    @error "workspace not initialized, please run setup(problemDict,bbSettingsDict,ssSettingsDict)"
  end
  return remove_constraints!(globalWorkspace,indices,
                                    suppressWarnings=suppressWarnings,
                                    suppressUpdate=suppressUpdate,
                                    localOnly=localOnly)
end


function permute_constraints_b(permutation::Array{Int,1},
                               suppressWarnings::Bool=false,
                               suppressUpdate::Bool=false,
                               localOnly::Bool=false)::Nothing
  if globalWorkspace isa NullWorkspace
    @error "workspace not initialized, please run setup(problemDict,bbSettingsDict,ssSettingsDict)"
  end

  return permute_constraints!(globalWorkspace,permutation,
                                     suppressWarnings=suppressWarnings,
                                     suppressUpdate=suppressUpdate,
                                     localOnly=localOnly)
end

function update_bounds_b(boundsDict::Dict,
                         suppressWarnings::Bool=false,
                         suppressUpdate::Bool=false,
                         localOnly::Bool=false)::Nothing
  if globalWorkspace isa NullWorkspace
    @error "workspace not initialized, please run setup(problemDict,bbSettingsDict,ssSettingsDict)"
  end

  boundsIn = Dict{Symbol,Array{Float64,1}}()
  for key in keys(boundsDict)
    boundsIn[Symbol(key)] = boundsDict[key]
  end
  return update_bounds!(globalWorkspace; boundsIn...,
                               suppressWarnings=suppressWarnings,
                               suppressUpdate=suppressUpdate,
                               localOnly=localOnly)
end


function append_problem_b(problemDict::Dict,
                          suppressWarnings::Bool=false,
                          suppressUpdate::Bool=false,
                          localOnly::Bool=false)::Nothing
  if globalWorkspace isa NullWorkspace
    @error "workspace not initialized, please run setup(problemDict,bbSettingsDict,ssSettingsDict)"
  end

  problemIn = Problem(problemDict)

  return append_problem!(globalWorkspace,problemIn,
                                suppressWarnings=suppressWarnings,
                                suppressUpdate=suppressUpdate,
                                localOnly=localOnly)
end


function integralize_variables_b(newDscIndices::Array{Int,1},
                                 newSos1Groups::Array{Int,1}=Int[],
                                 suppressWarnings::Bool=false,
                                 suppressUpdate::Bool=false,
                                 localOnly=false)::Nothing
  if globalWorkspace isa NullWorkspace
    @error "workspace not initialized, please run setup(problemDict,bbSettingsDict,ssSettingsDict)"
  end
  return integralize_variables!(globalWorkspace,newDscIndices,
                                suppressWarnings=suppressWarnings,
                                suppressUpdate=suppressUpdate,
                                localOnly=localOnly)
end

function update_objectiveCutoff_b(newCutoff::Float64,
                                  suppressWarnings::Bool=false,
                                  suppressUpdate::Bool=false,
                                  localOnly::Bool=false)::Nothing
  if globalWorkspace isa NullWorkspace
    @error "workspace not initialized, please run setup(problemDict,bbSettingsDict,ssSettingsDict)"
  end
  return update_objectiveCutoff!(globalWorkspace,newCutoff,
                                 suppressWarnings=suppressWarnings,
                                 suppressUpdate=suppressUpdate,
                                 localOnly=localOnly)
end
