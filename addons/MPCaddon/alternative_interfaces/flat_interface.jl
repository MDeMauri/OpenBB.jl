# @Author: Massimo De Mauri <massimo>
# @Date:   2019-06-19T16:34:10+02:00
# @Email:  massimo.demauri@gmail.com
# @Filename: ctypes_interface.jl
# @Last modified by:   massimo
# @Last modified time: 2021-01-27T00:34:03+01:00
# @License: LGPL-3.0
# @Copyright: {{copyright}}


######################## Shifts ########################
##### Common

### Shifts
# allowed modes are
#	"fullRTI" : problem shift and complete tree shift
#	"constraintsOnly": problem shift only (BB tree restarted from scratch)
# 	"warmStart": problem shift and partial tree shift used as heuristics
function mpc_shiftB(workspace::SupersolverWorkspace,shiftSteps::Int,newTailCnssDict::Dict,newTailCostsDict::Dict,
					referenceSolution::Array{Float64,1}=Float64[],measuredState::Array{Float64,1}=Float64[],
					mode::String="fullRTI",subsolverIterations::Int=0,
					suppressErrors::Bool=false)::Nothing

		return mpc_shift!(workspace,shiftSteps,ObjectiveFunction(newTailCostsDict),ConstraintSet(newTailCnssDict),
						  referenceSolution=referenceSolution,measuredState=measuredState,
						  mode=Symbol(mode),subsolverIterations=subsolverIterations,
						  suppressErrors=suppressErrors)
end

### workspace manipulation
function insert_constraints_ocB(workspace::SupersolverWorkspace,constraintsDict::Dict;suppressErrors::Bool=false)::Nothing
		return insert_constraints_oc!(workspace,ConstraintSet(constraintsDict),suppressErrors=suppressErrors)
end


##### BB Only

##### HBB Only
function mpc_shiftB(workspace::HBBworkspace,shiftSteps::Int,shiftedProblemDict::Dict,
					referenceSolution::Array{Float64,1}=Float64[],measuredState::Array{Float64,1}=Float64[],
					mode::String="fullRTI",subsolverIterations::Int=0,suppressErrors::Bool=false)::Nothing

		return mpc_shift!(workspace,shiftSteps,Problem(shiftedProblemDict),
						  referenceSolution=referenceSolution,measuredState=measuredState,
						  mode=Symbol(mode),subsolverIterations=subsolverIterations,suppressErrors=suppressErrors)
end
