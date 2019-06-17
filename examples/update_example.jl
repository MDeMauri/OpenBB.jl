# @Author: Massimo De Mauri <massimodemauri>
# @Date:   2019-02-06T16:22:14+01:00
# @Email:  massimo.demauri@gmail.com
# @Project: OpenBB
# @Filename:
# @Last modified by:   massimo
# @Last modified time: 2019-05-06T18:49:24+02:00
# @License: LGPL-3.0
# @Copyright: {{copyright}}

using OpenBB
using SparseArrays
using LinearAlgebra

clearconsole()

# subsolver = "osqp"
subsolverSettings = OpenBB.OSQPsettings()




println("---------------------------------------------------------------------------------------------------------------------")
println("First run")
println("---------------------------------------------------------------------------------------------------------------------")
# Basic usage of OpenBB for mixed-integer quadratic problems
problem = OpenBB.Problem(objFun=OpenBB.QuadraticObjective(Q=Matrix(2.0I,4,4,),L=[-.5,0.,0.,0.]),
                         cnsSet=OpenBB.LinearConstraintSet(A=ones(0,4),loBs=Float64[],upBs=Float64[]),
                         varSet=OpenBB.VariableSet(loBs=[-5.;-Infs(3)],upBs=[ 5.;Infs(3)],val=zeros(4),dscIndices=[1]))


workspace = OpenBB.setup(OpenBB.BBsettings(verbose=true,dynamicMode=true),subsolverSettings)
OpenBB.append_problem!(workspace,problem)

result = OpenBB.solve!(workspace)


println("")
println("---------------------------------------------------------------------------------------------------------------------")
println("Add constraints")
println("---------------------------------------------------------------------------------------------------------------------")
# add some linear contraints
OpenBB.append_constraints!(workspace,ones(1,4),[1.],[1.])
result2 = OpenBB.solve!(workspace)


println("")
println("---------------------------------------------------------------------------------------------------------------------")
println("Append problem")
println("---------------------------------------------------------------------------------------------------------------------")
# Basic usage of OpenBB for mixed-integer quadratic problems
problem2 = OpenBB.Problem(objFun=OpenBB.QuadraticObjective(Q=sparse([1,2,3,4],[1,2,3,4],[2.,4.,6.,8.]),L=[2.,2.,2.,2.]),
                         cnsSet=OpenBB.LinearConstraintSet(A=ones(1,4),loBs=[1.],upBs=[1.]),
                         varSet=OpenBB.VariableSet(loBs=[-5.;-Infs(3)],upBs=[ 5.;Infs(3)],val=zeros(4),dscIndices=[1]))

OpenBB.append_problem!(workspace,problem2)
result3 = OpenBB.solve!(workspace)


println("")
println("---------------------------------------------------------------------------------------------------------------------")
println("Update bounds")
println("---------------------------------------------------------------------------------------------------------------------")
OpenBB.update_bounds!(workspace;varLoBs=[1.,0.,0.,0.,1.,0.,0.,0.])
result4 = OpenBB.solve!(workspace)
