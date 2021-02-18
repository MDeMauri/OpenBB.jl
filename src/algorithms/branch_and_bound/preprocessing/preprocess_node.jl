# @Author: Massimo De Mauri <massimo>
# @Date:   2020-11-06T18:02:58+01:00
# @Email:  massimo.demauri@gmail.com
# @Filename: preprocess_node.jl
# @Last modified by:   massimo
# @Last modified time: 2021-01-25T13:48:42+01:00
# @License: LGPL-3.0
# @Copyright: {{copyright}}


function preprocess_node!(node::BBnode,workspace::BBworkspace,varsToCheck::Vector{Int64};withBoundsPropagation::Bool=true)::Bool

       # prepare output
       newUpdatedVars::Array{Int64, 1} = []
       feasible = true

      # Check bounds propagation
      if feasible && withBoundsPropagation

          # if the constraint set is non-linear (convex) obtain a linear relaxation. Do nothing otherwise.
          linearRlx = linearRelaxation(workspace.problem.cnsSet,node.primal)

          feasible, updatedVars = bounds_propagation!(node,linearRlx.A,
                                                           workspace.problem.varSet.dscIndices,
                                                           varsToCheck)
          union!(newUpdatedVars,collect(updatedVars))


          # Check SOS1 constraints
          if feasible && withBoundsPropagation
              feasible, updatedVars = preprocess_sos1!(newUpdatedVars,
                                                       get_sos1Groups(workspace),
                                                       workspace.problem.varSet.dscIndices,
                                                       node.varLoBs, node.varUpBs)
              union!(newUpdatedVars, collect(updatedVars))
          end

           # Recursively apply the preprocessing to the new updated variables
          if feasible && length(newUpdatedVars) > 0
              feasible = preprocess_node!(node,workspace,newUpdatedVars,withBoundsPropagation=withBoundsPropagation)
          end
      end

   return feasible
end


function preprocess_rows!(node::BBnode, workspace::BBworkspace, rowsToCheck::Vector{Int64})::Bool
      feasible, updatedVars = bounds_propagation!(
          Set(rowsToCheck),
          linearRelaxation.A,
          node.cnsLoBs, node.cnsUpBs,
          node.varLoBs, node.varUpBs,
          workspace.problem.varSet.dscIndices
      )
      updatedVarsArray = collect(updatedVars)::Array{Int64, 1}
      return preprocess_node!(node, workspace, updatedVarsArray, withBoundsPropagation=false)
end
