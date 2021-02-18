# @Author: Massimo De Mauri <massimo>
# @Date:   2020-11-24T17:34:57+01:00
# @Email:  massimo.demauri@gmail.com
# @Filename: _include.jl
# @Last modified by:   massimo
# @Last modified time: 2020-12-27T13:56:02+01:00
# @License: LGPL-3.0
# @Copyright: {{copyright}}


# abstract type
abstract type HBBmipStepWorkspace <: AbstractWorkspace end
abstract type HBBmipStepSettings <: AbstractSettings end


# these are generic fuctions that can be overloaded to have a specific implementation
# for a specific type of MIP step
function append_constraints!(workspace::HBBmipStepWorkspace,constraintSet::ConstraintSet;suppressErrors::Bool=false)::Nothing
    return append_constraints!(workspace.mipSolverWS,constraintSet,suppressErrors=suppressErrors)
end

function insert_constraints!(workspace::HBBmipStepWorkspace,constraintSet::ConstraintSet;suppressErrors::Bool=false)::Nothing
    return insert_constraints!(workspace.mipSolverWS,constraintSet,suppressErrors=suppressErrors)
end

function remove_constraints!(workspace::HBBmipStepWorkspace,indices::Union{Vector{Int},UnitRange{Int}};suppressErrors::Bool=false)::Nothing
    return remove_constraints!(workspace.mipSolverWS,indices,suppressErrors=suppressErrors)
end

function update_bounds!(workspace::HBBmipStepWorkspace;
						varLoBs::Vector{Float}=Vector{Float}(),
						varUpBs::Vector{Float}=Vector{Float}(),
                        cnsLoBs::Vector{Float}=Vector{Float}(),
                        cnsUpBs::Vector{Float}=Vector{Float}(),
                        suppressErrors::Bool=false)::Nothing

	return update_bounds!(workspace.mipSolverWS,varLoBs,varUpBs,cnsLoBs,cnsUpBs,suppressErrors=suppressErrors)
end

function append_problem!(workspace::HBBmipStepWorkspace,problem::Problem;suppressErrors::Bool=false)::Nothing
	return append_problem!(workspace.mipSolverWS,problem,suppressErrors=suppressErrors)
end


#
function get_variableBounds(workspace::HBBmipStepWorkspace)::Tuple{Vector{Float},Vector{Float}}
    return get_bounds(workspace.mipSolverWS)
end

#
function get_constraintBounds(workspace::HBBmipStepWorkspace)::Tuple{Vector{Float},Vector{Float}}
    return get_bounds(workspace.mipSolverWS)
end


function update_objectiveCutoff!(workspace::HBBmipStepWorkspace,newCutoff::Float)::Nothing
    if newCutoff < workspace.relaxationCutoff
        update_objectiveCutoff!(workspace.mipSolverWS,newCutoff)
        workspace.relaxationCutoff = newCutoff
    end
    return
end

function append_relaxations!(workspace::HBBmipStepWorkspace,newRelaxations::LinearConstraintSet)::Nothing
    append_constraints!(workspace.mipSolverWS,newRelaxations)
    return
end

function addto_blackList!(workspace::HBBmipStepWorkspace,assignment::Vector{Float})::Nothing
	addto_blackList!(workspace.mipSolverWS,assignment)
	return
end



include("./OAmipStep.jl")
