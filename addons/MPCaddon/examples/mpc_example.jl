# @Author: Massimo De Mauri <massimo>
# @Date:   2019-03-05T20:29:28+01:00
# @Email:  massimo.demauri@gmail.com
# @Filename: exampleMPC.jl
# @Last modified by:   massimo
# @Last modified time: 2021-02-14T19:45:05+01:00
# @License: LGPL-3.0
# @Copyright: {{copyright}}
using Revise
using Distributed
using SparseArrays
using LinearAlgebra
using OpenBB

# construct the short time horizon problem
horizonLen = 25
initialState = 10.5
mpcSteps = 10
numStates = 1
numControls = 3
numVars = numStates + numControls

# objective at each step
Q_ = sparse(Array{Float64,2}(2*I,numVars,numVars))
Q_[3,3] = 0.

# build full objective matrices
Q = spzeros(0,0)
for s in 1:horizonLen
    global Q = vcat(hcat(Q,spzeros(numVars*(s-1),numVars)),hcat(spzeros(numVars,numVars*(s-1)),Q_))
end
Q = vcat(hcat(Q,spzeros(size(Q,1),1)),hcat(spzeros(1,size(Q,1)),0.0))
L = zeros(numVars*horizonLen+1)

# constraints at each step
pathCns = [0.,1.,1.,1.,0.]'
dynamicCns = [-1.,1.,0.,-1,1.]'
A_ = sparse(vcat(pathCns,dynamicCns))

# build full constraint set
A = hcat([1],zeros(1,numVars*horizonLen)) #imposes the initial state
cnsLoBs = [initialState]
cnsUpBs = [initialState]
for s in 1:horizonLen
    global A = vcat(A,hcat(zeros(size(A_,1),(s-1)*numVars),A_,zeros(size(A_,1),(horizonLen-s)*numVars)))
    global cnsLoBs = vcat(cnsLoBs,[1.,0.])
    global cnsUpBs = vcat(cnsUpBs,[1.,0.])
end


varLoBs = vcat(repeat([-100.,0.,0.,0.],horizonLen),[-100.])
varUpBs = vcat(repeat([ 100.,1.,1., 1.],horizonLen),[100.])
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
workspace = OpenBB.setup(problem,OpenBB.BBsettings(verbose=true,
                                                   optimalControlInfo=(numStates,0,numControls),
                                                   conservativismLevel=2,
                                                   withBoundsPropagation=true,
                                                   numProcesses=1,
                                                   subsolverSettings=OpenBB.OSQPsettings()))

# define the tail cost and constraints for the problem
tailCnsSet = OpenBB.LinearConstraintSet(A=A[end-1:end,end-2*numStates-numControls+1:end],loBs=cnsLoBs[end-1:end],upBs=cnsUpBs[end-1:end])
tailObjFun = OpenBB.QuadraticObjective(Q=Q[end-2*numStates-numControls+1:end,end-2*numStates-numControls+1:end],L=L[end-2*numStates-numControls+1:end])

starting_time = time()
for s in 1:mpcSteps
    println("\n=======================================================",
              "=======================================================")
    println("Iteration ",s,":")

    # solve the mpc subproblem
    OpenBB.solve!(workspace)
    solution = OpenBB.get_best_feasible_node(workspace)

    # print results
    println()
    println("Results:")
    println("   - cost: ",solution.objUpB)

    stateTrajectory = Array{Float64,2}(undef,0,horizonLen+1)
    for k in 1:numStates
        stateTrajectory = vcat(stateTrajectory,transpose(solution.primal[k:numVars:end]))
    end
    println("   - state trajectory:")
    for k in 1:size(stateTrajectory,1)
        println("      x",k," = ",@. round(stateTrajectory[k,:],digits=3))
    end

    controlAction = Array{Float64,2}(undef,0,horizonLen)
    for k in numStates+1:numVars
        controlAction = vcat(controlAction,transpose(solution.primal[k:numVars:end]))
    end
    println("   - control action: ")
    for k in 1:size(controlAction,1)
        println("      c",k," = ",@. round(controlAction[k,:],digits=3))
    end

    # shift the problem one step forward
    if s < mpcSteps
        OpenBB.mpc_shift!(workspace,1,tailObjFun,tailCnsSet,mode=:fullRTI,referenceSolution=solution.primal,approximateBellmanConstant=0.0)
    end
end
println("Total time: ",time()-starting_time)
