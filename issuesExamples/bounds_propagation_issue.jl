# @Author: Massimo De Mauri <massimo>
# @Date:   2019-09-02T14:19:22+02:00
# @Email:  massimo.demauri@gmail.com
# @Filename: bounds_propagation_issue.jl
# @Last modified by:   massimo
# @Last modified time: 2019-09-06T18:22:12+02:00
# @License: LGPL-3.0
# @Copyright: {{copyright}}


using Distributed
using SparseArrays
using LinearAlgebra
using OpenBB


# construct the short time horizon problem
horizonLen = 30
initialState = 0.5
mpcSteps = 3
numStates = 1
numControls = 3
numVars = numStates + numControls

# objective
Q_ = sparse(Array{Float64,2}(2*I,numVars,numVars))
Q_[3,3] = 0.
Q = spzeros(0,0)
for s in 1:horizonLen
    global Q = vcat(hcat(Q,spzeros(numVars*(s-1),numVars)),hcat(spzeros(numVars,numVars*(s-1)),Q_))
end
Q = vcat(hcat(Q,spzeros(size(Q,1),1)),hcat(spzeros(1,size(Q,1)),0.5))
L = zeros(numVars*horizonLen+1)

# constraints
A_ = transpose([[1.,1.,0.,-1,-1.] [0.,1.,1.,1.,0.]])
A = sparse(A_)
for s in 2:horizonLen
    global A = vcat(hcat(A,spzeros(size(A_,1)*(s-1),numVars)),hcat(spzeros(size(A_,1),numVars*(s-1)),A_))
end
cnsUpBs = cnsLoBs = repeat([0.,1.],horizonLen)

varLoBs = vcat(initialState,repeat([0.,0.,0.,-100.],horizonLen))
varUpBs = vcat(initialState,repeat([1.,1., 1., 100.],horizonLen))
vars_val = vcat(initialState,repeat([0.,1.,0.,initialState],horizonLen))
dscIndices = [numVars*(k-1) + i + numStates for k in 1:horizonLen for i in 1:numControls]
sos1Groups = [k for k in 1:horizonLen for i in 1:numControls]
# sos1Groups = [-1 for k in 1:horizonLen for i in 1:numControls]


# build problem componets
cnss = OpenBB.LinearConstraintSet(A=A,loBs=cnsLoBs,upBs=cnsUpBs)
objf = OpenBB.QuadraticObjective(Q=Q,L=L)
vars = OpenBB.VariableSet(loBs=varLoBs,upBs=varUpBs,vals=vars_val,dscIndices=dscIndices,sos1Groups=sos1Groups)


# build the problem
problem = OpenBB.Problem(objFun=objf,cnsSet=cnss,varSet=vars)
workspace = OpenBB.setup(problem,
                         OpenBB.BBsettings(verbose=true,
                                           interactiveMode=true,
                                           expansionPriorityRule=(OpenBB.lower_pseudoObj,),
                                           branchingPriorityRule=(OpenBB.pseudoIncrements_geomean,),
                                           pseudoCostsInitialization=(OpenBB.initialize_to_constant!,1e-4),
                                           withBoundsPropagation=false,
                                           numProcesses=1),
                         OpenBB.GUROBIsettings())

workspacePreprocess = OpenBB.setup(problem,
                         OpenBB.BBsettings(verbose=true,
                                           interactiveMode=true,
                                           expansionPriorityRule=(OpenBB.lower_pseudoObj,),
                                           branchingPriorityRule=(OpenBB.pseudoIncrements_geomean,),
                                           pseudoCostsInitialization=(OpenBB.initialize_to_constant!,1e-4),
                                           withBoundsPropagation=true,
                                           numProcesses=1),
                         OpenBB.GUROBIsettings())

# solve
OpenBB.solve!(workspacePreprocess)
OpenBB.solve!(workspace)
