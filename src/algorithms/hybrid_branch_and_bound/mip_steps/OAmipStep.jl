# @Author: Massimo De Mauri <massimo>
# @Date:   2020-11-24T17:35:43+01:00
# @Email:  massimo.demauri@gmail.com
# @Filename: OAmipStep.jl
# @Last modified by:   massimo
# @Last modified time: 2021-02-16T14:19:47+01:00
# @License: LGPL-3.0
# @Copyright: {{copyright}}

mutable struct OAmipStepSettings <: HBBmipStepSettings
    mipSettings::BBsettings
end

# named constructor
function OAmipStepSettings(;mipSettings::BBsettings=BBsettings())::OAmipStepSettings
	return OAmipStepSettings(mipSettings)
end


mutable struct OAmipStepWorkspace <: HBBmipStepWorkspace
    # info for the nlpSubsolver
    mipSolverWS::BBworkspace
	# map the relaxations into the original constraints
	relaxationsIDs::Vector{Int}
    # info for the step
    relaxationCutoff::Float
    # tells if the workspace needs an update
    outdated::Bool
end



function setup(relaxation::Problem,relaxationsIDs::Vector{Int},stepSettings::OAmipStepSettings)::OAmipStepWorkspace

    @assert relaxation.cnsSet isa LinearConstraintSet &&
            (relaxation.objFun isa QuadraticObjective || relaxation.objFun isa LinearObjective)

    # just create a branch and bound process for the given relaxation
    mipSolverWS = setup(relaxation,stepSettings.mipSettings)

    return OAmipStepWorkspace(mipSolverWS,relaxationsIDs,Inf,false)
end

function make_outdated!(workspace::OAmipStepWorkspace)::Nothing
    workspace.outdated = true
    make_outdated!(workspace.mipSolverWS)
    return
end

function update!(workspace::OAmipStepWorkspace)::Nothing
    update!(workspace.mipSolverWS)
    return
end


function solve!(workspace::OAmipStepWorkspace)::Tuple{String,Vector{BBnode},Float}

    # perform branch and bound
    solve!(workspace.mipSolverWS)

    # get status and objective value
    if workspace.mipSolverWS.status.description in ["optimalSolutionFound","interrupted"]
        solutions = get_all_feasible_nodes(workspace.mipSolverWS)
        if workspace.mipSolverWS.sharedMemory isa NullSharedMemory
            bestBound = workspace.mipSolverWS.status.objLoB
        else
            bestBound = minimum(workspace.mipSolverWS.sharedMemory.localObjLoBs)
        end
    elseif workspace.mipSolverWS.status.description in ["cutoffInfeasible","blackListInfeasible","infeasible"]
        solutions = BBnode[]
        bestBound = Inf
    elseif workspace.mipSolverWS.status.description == "unreliableSolution"
        solutions = get_all_feasible_nodes(workspace.mipSolverWS)
        if workspace.mipSolverWS.sharedMemory isa NullSharedMemory
            bestBound = workspace.mipSolverWS.status.objLoB
        else
            bestBound = minimum(workspace.mipSolverWS.sharedMemory.localObjLoBs)
        end
        error("handling of unreliable solutions still to implement")
    end

    return workspace.mipSolverWS.status.description,solutions,bestBound

end

## inspect functions
function get_mipSettings(workspace::OAmipStepWorkspace)::BBsettings
	return workspace.mipSolverWS.settings
end

function get_relaxationsSet(workspace::OAmipStepWorkspace)::LinearConstraintSet
	return workspace.mipSolverWS.problem.cnsSet
end

function get_relaxationsIDs(workspace::OAmipStepWorkspace)::Vector{Int}
	return workspace.relaxationsIDs
end

## update functions

function rename_relaxations!(workspace::OAmipStepWorkspace,newIDs::Vector{Int})::Nothing
	workspace.relaxationsIDs = newIDs[workspace.relaxationsIDs]
	return
end

# these are generic fuctions that can be overloaded to have a specific implementation
# for a specific type of MIP step\
function append_relaxations!(workspace::OAmipStepWorkspace,cnsSet::ConstraintSet,IDs::Vector{Int};suppressErrors::Bool=false)::Nothing
	append!(workspace.relaxationsIDs,IDs)
    return append_constraints!(workspace.mipSolverWS,cnsSet,suppressErrors=suppressErrors)
end


function remove_constraints!(workspace::OAmipStepWorkspace,IDs::Union{Vector{Int},UnitRange{Int}};suppressErrors::Bool=false)::Nothing
	indices = findall(x->x in IDs,workspace.relaxationsIDs)
    remove_constraints!(workspace.mipSolverWS,indices,suppressErrors=suppressErrors)
	deleteat!(workspace.relaxationsIDs,indices)
	return
end

function update_bounds_variables!(workspace::OAmipStepWorkspace;
								  loBs::Vector{Float}=Vector{Float}(),
								  upBs::Vector{Float}=Vector{Float}(),
                        		  suppressErrors::Bool=false)::Nothing

	return update_bounds!(workspace.mipSolverWS,varLoBs=loBs,varUpBs=upBs,suppressErrors=suppressErrors)
end

function update_bounds_constraints!(workspace::OAmipStepWorkspace;
									oldLoBs::Vector{Float}=Vector{Float}(),
									oldUpBs::Vector{Float}=Vector{Float}(),
                        			newLoBs::Vector{Float}=Vector{Float}(),
                        			newUpBs::Vector{Float}=Vector{Float}(),
                        			suppressErrors::Bool=false)::Nothing

	@assert length(oldLoBs) == length(newLoBs) && length(oldUpBs) == length(newUpBs)
	deltaLoBs = if !isempty(oldLoBs) newLoBs - oldLoBs else Float[] end
	deltaUpBs = if !isempty(oldUpBs) newUpBs - oldUpBs else Float[] end
	rlxLoBs = if !isempty(deltaLoBs) workspace.mipSolverWS.problem.cnsSet.loBs + deltaLoBs[workspace.relaxationsIDs] else Float[] end
	rlxUpBs = if !isempty(deltaUpBs) workspace.mipSolverWS.problem.cnsSet.upBs + deltaUpBs[workspace.relaxationsIDs] else Float[] end
	return update_bounds!(workspace.mipSolverWS,cnsLoBs=rlxLoBs,cnsUpBs=rlxUpBs,suppressErrors=suppressErrors)
end


function append_problem!(workspace::OAmipStepWorkspace,problem::Problem,newRelaxationsIDs::Vector{Int};suppressErrors::Bool=false)::Nothing
	append!(workspace.linearizationsIDs,newRelaxationsIDs)
	return append_problem!(workspace.mipSolverWS,problem,suppressErrors=suppressErrors)
end


function fix_variables!(workspace::OAmipStepWorkspace,indices::Union{Vector{Int},UnitRange{Int}},values::Vector{Float};
						removeFixedVariables::Bool=false,suppressErrors::Bool=false)::Nothing
	return fix_variables!(workspace.mipSolverWS,indices,values,removeFixedVariables=removeFixedVariables,suppressErrors=suppressErrors)
end

function integralize_variables!(workspace::OAmipStepWorkspace,indices::Vector{Int})::Nothing
	return integralize_variables!(workspace.mipSolverWS)
end

#
function get_variableBounds(workspace::OAmipStepWorkspace)::Tuple{Vector{Float},Vector{Float}}
    return get_variableBounds(workspace.mipSolverWS)
end

#
function get_relaxationsBounds(workspace::OAmipStepWorkspace)::Tuple{Vector{Float},Vector{Float}}
    return get_constraintBounds(workspace.mipSolverWS)
end


function update_objectiveCutoff!(workspace::OAmipStepWorkspace,newCutoff::Float)::Nothing
    if newCutoff < workspace.relaxationCutoff
        update_objectiveCutoff!(workspace.mipSolverWS,newCutoff)
        workspace.relaxationCutoff = newCutoff
    end
    return
end

function permute_constraints!(workspace::OAmipStepWorkspace,permutation::Vector{Int})::Nothing
	inversePermutation = invperm(permutation)
	workspace.relaxationsIDs .= inversePermutation[workspace.relaxationsIDs]
	return
end

function permute_relaxations!(workspace::OAmipStepWorkspace,permutation::Vector{Int})::Nothing
	@assert length(permutation) == length(workspace.relaxationsIDs)
	permute!(workspace.relaxationsIDs,permutation)
	return permute_constraints!(workspace.mipSolverWS,permutation)
end



function addto_blackList!(workspace::OAmipStepWorkspace,assignment::Vector{Float})::Nothing
	addto_blackList!(workspace.mipSolverWS,assignment)
	return
end

function get_blacklist(workspace::OAmipStepWorkspace)::BBblackList
	return workspace.mipSolverWS.blacklist
end

function clear_blackList!(workspace::OAmipStepWorkspace)::Nothing
	return clear_blackList!(workspace.mipSolverWS)
end

function Base.filter!(selector::Function,workspace::OAmipStepWorkspace,suppressErrors::Bool)::Nothing
	return filter!(selector,workspace.mipSolverWS,suppressErrors=suppressErrors)
end
