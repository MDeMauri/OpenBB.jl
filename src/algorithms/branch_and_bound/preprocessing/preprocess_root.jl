# @Author: Massimo De Mauri <massimo>
# @Date:   2020-11-06T17:52:26+01:00
# @Email:  massimo.demauri@gmail.com
# @Filename: root_preprocessing.jl
# @Last modified by:   massimo
# @Last modified time: 2021-02-12T18:05:47+01:00
# @License: LGPL-3.0
# @Copyright: {{copyright}}


function preprocess_root!(node::BBnode,workspace::BBworkspace)::Bool
   # all nodes are presumed feasible until until proved infeasible
   feasible = true
   if workspace.problem.cnsSet isa LinearConstraintSet
      feasible = preprocess_gcd!(workspace.problem.cnsSet.A,
                                 node.cnsLoBs, node.cnsUpBs,
                                 node.varLoBs, node.varUpBs,
                                 workspace.problem.varSet.dscIndices)
   end
   return feasible
end
