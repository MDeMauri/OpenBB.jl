# @Author: Massimo De Mauri <massimodemauri>
# @Date:   2019-02-06T16:22:14+01:00
# @Email:  massimo.demauri@gmail.com
# @Project: OpenBB
# @Filename:
# @Last modified by:   massimo
# @Last modified time: 2021-02-12T13:11:29+01:00
# @License: LGPL-3.0
# @Copyright: {{copyright}}

using OpenBB
using SparseArrays
using LinearAlgebra

OpenBB.load_subsolver_interface("OSQP")
subsolverSettings = OpenBB.OSQPsettings()

# Basic usage of OpenBB for mixed-integer quadratic problems
problem = OpenBB.Problem(objFun=OpenBB.QuadraticObjective(Q=Matrix(2.0I,4,4,),L=[-.5,0.,0.,0.]),
                         cnsSet=OpenBB.LinearConstraintSet(A=ones(1,4),loBs=[1.],upBs=[1.]),
                         varSet=OpenBB.VariableSet(loBs=[-5.;-Infs(3)],upBs=[ 5.;Infs(3)],vals=zeros(4),dscIndices=[1]))

workspace = OpenBB.setup(problem,OpenBB.BBsettings(verbose=true,conservativismLevel=1,numProcesses=1,subsolverSettings=subsolverSettings))
result = OpenBB.solve!(workspace)
