# @Author: Massimo De Mauri <massimo>
# @Date:   2020-11-24T19:02:25+01:00
# @Email:  massimo.demauri@gmail.com
# @Filename: hybridBranchAndBound.jl
# @Last modified by:   massimo
# @Last modified time: 2021-02-12T14:29:14+01:00
# @License: LGPL-3.0
# @Copyright: {{copyright}}

using Revise
using SparseArrays
using LinearAlgebra
using OpenBB

OpenBB.load_subsolver_interface("IPOPT")
OpenBB.load_subsolver_interface("OSQP")

# create mixed integer non-linear problem (with linear/quadratic objective)
varSet = OpenBB.VariableSet(vals=zeros(4),loBs=-Infs(4),upBs=Infs(4),dscIndices=[1,2])
objFun = OpenBB.QuadraticObjective(Q=Matrix{Float64}(I,4,4),L=zeros(4))

function evalCnss(x)
    return [ x[1] + x[2],
             x[3]^2 + x[4]^2,
             x[3] - 10.0*x[1],
             x[4] -  5.0*x[2] ]
end
cnsLoBs = [1.0, -Inf, 0.0, 0.0]
cnsUpBs = [1.0,  1e2, Inf, Inf]

function evalJcb(x)
    return Matrix([    1.0         1.0     0.0         0.0         ;
                       0.0         0.0     2.0*x[3]    2.0*x[4]    ;
                       -10.0       0.0     1.0         0.0         ;
                       0.0         -5.0    0.0         1.0         ])
end
jcbSparsity = sparse(Matrix([   true      true    false   false   ;
                                false     false   true    true    ;
                                true      false   true    false   ;
                                false     true    false   true    ]))

evalHes = Vector{Function}(undef,4)
evalHes = x->[  zeros(4,4),
                zeros(4,4),
                vcat(zeros(2,4),hcat(zeros(2,2),Matrix(2.0I,2,2))),
                zeros(4,4)]

hesSparsity = [sparse(falses(4,4)),
              sparse(falses(4,4)),
              sparse(vcat(falses(2,4),hcat(falses(2,2),Matrix(true*I,2,2)))),
              sparse(falses(4,4))]

cnsSet = OpenBB.ConvexConstraintSet(evalVal=evalCnss,evalJcb=evalJcb,evalHes=evalHes,
                                       jcbSparsity=jcbSparsity,hesSparsity=hesSparsity,
                                       typeJcb=Matrix{Float64},typeHes=Matrix{Float64},
                                       loBs=cnsLoBs,upBs=cnsUpBs)

problem = OpenBB.Problem(varSet=varSet,cnsSet=cnsSet,objFun=objFun)

workspace = OpenBB.setup(problem,OpenBB.HBBsettings(verbose=true,nlpSettings=OpenBB.IPOPTsettings(),mipSettings=OpenBB.BBsettings(verbose=true,subsolverSettings=OpenBB.OSQPsettings())))
OpenBB.solve!(workspace)

# get solution
optimalNode = OpenBB.get_best_feasible_node(workspace)
@info optimalObjective = optimalNode.objUpB
@info optimalSolution = optimalNode.primal
