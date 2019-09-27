# @Author: Wim Van Roy
# @Date:   2019-03-11T12:27:34+01:00
# @Email:  wim.vanroy@kuleuven.be
# @Filename: preprocessing.jl
# @Last modified by:   massimo
# @Last modified time: 2019-09-26T12:51:06+02:00
# @License: LGPL-3.0
# @Copyright: {{copyright}}

struct InfeasibleError <: Exception
   msg::String
end

include("./linear_bounds_propagation.jl")
include("./gcd.jl")
include("./sos.jl")


function preprocess!(node::BBnode, workspace::BBworkspace, varsToCheck::Array{Int64,1};
                     withBoundsPropagation::Bool=true)::Bool

   newUpdatedVars = Array{Int64, 1}

   feasible = true
   if withBoundsPropagation
      if 0 in varsToCheck
          feasible = preprocess_gcd!(
             get_linearConstraints(workspace.problem.cnsSet),
             node.cnsLoBs, node.cnsUpBs,
             node.varLoBs, node.varUpBs,
             workspace.problem.varSet.dscIndices,
          )
      end

      if feasible
          feasible, updatedVars = OpenBB.bounds_propagation!(
            node,
            get_linearConstraints(workspace.problem.cnsSet),
            workspace.problem.varSet.dscIndices,
            varsToCheck
          )
          newUpdatedVars = unique(vcat(newUpdatedVars, updatedVars))
      end

      if false
          # As we branch on sos1, it is not worth to do this!
          feasible, updatedVars = preprocess_sos1!(
            varsToCheck,
            get_sos1Groups(workspace),
            workspace.problem.varSet.dscIndices,
            node.varLoBs, node.varUpBs
          )
          newUpdatedVars = unique(vcat(newUpdatedVars, updatedVars))
      end

      # TO FIX
      # if feasible
      #     feasible = preprocess!(node, workspace, newUpdatedVars, withBoundsPropagation)
      # end
   end

   return feasible
end
