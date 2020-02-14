# @Author: Wim Van Roy <wimvr93>
# @Date:   Fr Feb 14 11:14:25 CET 2020
# @Email:  wim.vr@hotmail.com
# @Project: OpenBB
# @Filename:
# @Last modified by:   wimvr93
# @Last modified time: Fr Feb 14 11:14:58 CET 2020
# @License: LGPL-3.0
# @Copyright: {{copyright}}

using Revise
using MAT
using SparseArrays
using LinearAlgebra
using OpenBB

file = matopen(Base.source_dir()*"/problem_matrices.mat")

Q = read(file, "P")
L = convert(Array{Float64,1},read(file, "q")[:])
A = read(file, "A")
cnsLoBs = read(file, "lbg")[:]
cnsUpBs = read(file, "ubg")[:]
varLoBs = read(file, "lbw")[:]
varUpBs = read(file, "ubw")[:]
dscIndices = convert(Array{Int64, 1},read(file, "discrete_idx")[:]) # Discrete indices (0 or 1)

varVals = copy(varLoBs)

cnss = OpenBB.LinearConstraintSet(A=A,loBs=cnsLoBs,upBs=cnsUpBs)
# OpenBB.oc_sort!(cnss)

problem = OpenBB.Problem(objFun=OpenBB.QuadraticObjective(Q=sparse(zeros(length(L), length(L))),L=L),
                         cnsSet=cnss,
                         varSet=OpenBB.VariableSet(loBs=varLoBs,upBs=varUpBs,vals=varVals,dscIndices=dscIndices))

workspace = OpenBB.setup(
                    problem,
                    OpenBB.BBsettings(withBoundsPropagation=false,
                                      verbose=true,
                                      interactiveMode=true,
                                      numProcesses=1),
                    OpenBB.GUROBIsettings())
result = OpenBB.solve!(workspace)
