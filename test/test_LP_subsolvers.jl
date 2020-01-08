# @Author: Massimo De Mauri <massimo>
# @Date:   2020-01-08T16:52:49+01:00
# @Email:  massimo.demauri@gmail.com
# @Filename: test_LP_subsolvers.jl
# @Last modified by:   massimo
# @Last modified time: 2020-01-08T16:57:52+01:00
# @License: LGPL-3.0
# @Copyright: {{copyright}}

using OpenBB
using LinearAlgebra
using SparseArrays

function test_LP_subsolver(subsolver)

    if subsolver == "CLP"
        subsolverSettings = OpenBB.CLPsettings()
    else
        error("Subsolver Unknown "*subsolver)
    end

    print(" - ")
    # create first problem
    problem = OpenBB.Problem(objFun=OpenBB.LinearObjective(L=[-.5,0.,0.,0.]),
                             cnsSet=OpenBB.LinearConstraintSet(A=ones(1,4),loBs=[0.0],upBs=Float64[0.0]),
                             varSet=OpenBB.VariableSet(loBs=[-5.;-Infs(3)],upBs=[ 5.;Infs(3)],vals=zeros(4),dscIndices=[1]))
    workspace = OpenBB.setup(problem,OpenBB.BBsettings(interactiveMode=true,verbose=false,statusInfoPeriod=0.01,numProcesses=1),subsolverSettings)
    result0 = OpenBB.solve!(workspace)

    # add some linear contraints
    OpenBB.append_constraints!(workspace,OpenBB.LinearConstraintSet(ones(1,4),[1.],[1.]))
    result1 = OpenBB.solve!(workspace)

    # Basic usage of OpenBB for mixed-integer quadratic problems
    problem2 = OpenBB.Problem(objFun=OpenBB.LinearObjective(L=[2.,2.,2.,2.]),
                             cnsSet=OpenBB.LinearConstraintSet(A=ones(1,4),loBs=[1.],upBs=[1.]),
                             varSet=OpenBB.VariableSet(loBs=[-5.;-Infs(3)],upBs=[ 5.;Infs(3)],vals=zeros(4),dscIndices=[1]))

    OpenBB.append_problem!(workspace,problem2)
    OpenBB.permute_constraints!(workspace,reverse(collect(1:OpenBB.get_numConstraints(workspace))))
    result2 = OpenBB.solve!(workspace)

    OpenBB.update_bounds!(workspace;varLoBs=[1.,0.,0.,0.,1.,0.,0.,0.])
    result3 = OpenBB.solve!(workspace)
    println(subsolver,": setup + solve + update, ok")
end

for subsolver_info in OpenBB.get_available_subsolvers()
    if subsolver_info[2] == "LP"
        test_LP_subsolver(subsolver_info[1])
    end
end
