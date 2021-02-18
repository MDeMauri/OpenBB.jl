# @Author: Massimo De Mauri <massimo>
# @Date:   2020-11-12T19:59:30+01:00
# @Email:  massimo.demauri@gmail.com
# @Filename: test_NLP_subsolvers.jl
# @Last modified by:   massimo
# @Last modified time: 2021-02-12T14:33:44+01:00
# @License: LGPL-3.0
# @Copyright: {{copyright}}

using OpenBB
using LinearAlgebra
using SparseArrays

function test_NLP_subsolver(subsolver)

    # build default settings for the solver
    subsolverSettings = OpenBB.eval(Symbol(subsolver,"settings"))()

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


    # test non-linear branch&bound
    bbWorkspace = OpenBB.setup(problem,OpenBB.BBsettings(subsolverSettings=subsolverSettings))
    OpenBB.solve!(bbWorkspace)
    optimalNode = OpenBB.get_best_feasible_node(bbWorkspace)
    @assert abs(13.0 - optimalNode.objUpB) <= bbWorkspace.settings.primalTolerance

    # test hybrid-branch&bound
    hbbWorkspace = OpenBB.setup(problem,OpenBB.HBBsettings(nlpSettings=subsolverSettings,mipSettings=OpenBB.BBsettings(subsolverSettings=subsolverSettings)))
    OpenBB.solve!(hbbWorkspace)
    optimalNode = OpenBB.get_best_feasible_node(hbbWorkspace)
    @assert abs(13.0 - optimalNode.objUpB) <= hbbWorkspace.settings.primalTolerance

    return
end

for subsolver_info in OpenBB.get_available_subsolvers()
    # test branch and bound
    if "NLP" in subsolver_info
        OpenBB.load_subsolver_interface(subsolver_info[1])
        test_NLP_subsolver(subsolver_info[1])
    end
end
