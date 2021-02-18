# @Author: Massimo De Mauri <massimo>
# @Date:   2019-06-19T18:15:50+02:00
# @Email:  massimo.demauri@gmail.com
# @Filename: test_flat_interface.jl
# @Last modified by:   massimo
# @Last modified time: 2021-02-12T17:08:15+01:00
# @License: LGPL-3.0

using LinearAlgebra
using SparseArrays
using OpenBB

OpenBB.include_flat_interface()


# define QP problem
problemDict = Dict{String,Dict{String,Any}}("objFun"=>Dict(),"cnsSet"=>Dict(),"varSet"=>Dict())
problemDict["objFun"] = Dict{String,Any}("type"=>"Quadratic","Q"=>2.0*speye(5),"L"=>zeros(5))
problemDict["objFun"]["Q"][3,3] = 0
dropzeros!(problemDict["objFun"]["Q"])
problemDict["cnsSet"] = Dict{String,Any}("type"=>"Linear","A"=>vcat([1. 1. 0. -1. -1.],[0. 1. 1. 1. 0.]),"loBs"=>[0.,1.],"upBs"=>[0.,1.])
problemDict["varSet"] = Dict{String,Any}("vals"=> [.5,0,1,0,.5],"loBs"=>[.5,0,0,0,-10],"upBs"=>[.5,1,1,1,10],"dscIndices"=>[2,3,4],"sos1Groups"=>[1,1,1])
bbSettings = Dict{String,Any}("subsolverSettings"=>Dict{String,Any}("subsolverName"=>"OSQP"),"verbose"=>false,"numProcesses"=>1,"conservativismLevel"=>1,"withBoundsPropagation"=>true)


## BB

# solve MIQP
workspace = OpenBB.setup("BB",problemDict,bbSettings)
OpenBB.solveB(workspace)

solution = OpenBB.get_best_feasible_node_dict(workspace)
@assert 0.5 - workspace.settings.primalTolerance <= OpenBB.get_best_feasible_node_dict(workspace)["objUpB"] <= 0.5 + workspace.settings.primalTolerance
OpenBB.get_all_feasible_nodes_dict(workspace)
OpenBB.get_best_node_dict(workspace)

@assert OpenBB.get_numVariables(workspace) == 5
@assert OpenBB.get_numConstraints(workspace) == 2
@assert OpenBB.get_numDiscrete(workspace) == 3

@assert OpenBB.get_constraints_dependency(workspace) == [[1,2,4,5],[2,3,4]]
@assert OpenBB.get_constraint_dependency(workspace,1) == [1,2,4,5]
@assert OpenBB.get_constraint_dependency(workspace,2) == [2,3,4]
@assert OpenBB.get_objective_dependency(workspace) == [1, 2, 4, 5]

@assert OpenBB.get_variableBounds(workspace) == ([0.5, 0.0, 0.0, 0.0, -10.0], [0.5, 1.0, 1.0, 1.0, 10.0])
@assert OpenBB.get_constraintBounds(workspace) == ([0.0, 1.0], [0.0, 1.0])
@assert OpenBB.get_status_dict(workspace)["description"] == "optimalSolutionFound"

OpenBB.resetB(workspace)
@assert OpenBB.get_numActiveNodes(workspace) == 1
@assert OpenBB.get_status_dict(workspace)["objUpB"] == Inf

OpenBB.clearB(workspace)
@assert OpenBB.get_numActiveNodes(workspace) == 0

OpenBB.resetB(workspace)
@assert OpenBB.get_numFeasibleNodes(workspace) == 0
@assert OpenBB.get_numActiveNodes(workspace) == 1

OpenBB.solveB(workspace)
OpenBB.append_constraintsB(workspace,problemDict["cnsSet"],true)
OpenBB.remove_constraintsB(workspace,[3,4],true)
OpenBB.permute_constraintsB(workspace,[2,1],true)
@assert OpenBB.get_constraintBounds(workspace) == ([1.0, 0.0], [1.0, 0.0])
OpenBB.update_boundsB(workspace,Dict("cnsLoBs"=>[0.0,0.0],"varLoBs"=>[0.5, 0.0, 0.0, 0.0, -100.0]),true)
@assert OpenBB.get_constraintBounds(workspace) == ([0.0, 0.0], [1.0, 0.0])
@assert OpenBB.get_variableBounds(workspace) == ([0.5, 0.0, 0.0, 0.0, -100.0], [0.5, 1.0, 1.0, 1.0, 10.0])
OpenBB.update_boundsB(workspace,Dict("cnsLoBs"=>[1.0,0.0]),true)


OpenBB.append_problemB(workspace,problemDict,true)
OpenBB.updateB(workspace)
OpenBB.solveB(workspace)
@assert 1.0 - workspace.settings.primalTolerance <= OpenBB.get_best_feasible_node_dict(workspace)["objUpB"] <= 1.0 + workspace.settings.primalTolerance


OpenBB.integralize_variablesB(workspace,[5])
OpenBB.solveB(workspace)
@assert OpenBB.get_status_dict(workspace)["description"] == "infeasible"



## HBB

# define NLP problem
problemDict = Dict{String,Dict{String,Any}}("objFun"=>Dict(),"cnsSet"=>Dict(),"varSet"=>Dict())
problemDict["objFun"] = Dict{String,Any}("type"=>"Quadratic","Q"=>2.0*speye(5),"L"=>zeros(5))
problemDict["objFun"]["Q"][3,3] = 0
dropzeros!(problemDict["objFun"]["Q"])
A = vcat([1. 1. 0. -1. -1.],[0. 1. 1. 1. 0.])
problemDict["cnsSet"] = Dict{String,Any}("type"=>"Convex",
                                         "evalVal"=>x->A*x,"evalJcb"=>x->A,"evalHes"=>x->[zeros(5,5),zeros(5,5)],
                                         "loBs"=>[0.,1.],"upBs"=>[0.,1.],"typeJcb"=>Matrix{Float64},"typeHes"=>Matrix{Float64},
                                         "jcbSparsity"=>OpenBB.sparsityMatrix(A),"hesSparsity"=>[spzeros(Bool,5,5),spzeros(Bool,5,5)])
problemDict["varSet"] = Dict{String,Any}("vals"=> [.5,0,1,0,.5],"loBs"=>[.5,0,0,0,-10],"upBs"=>[.5,1,1,1,10],"dscIndices"=>[2,3,4],"sos1Groups"=>[1,1,1])

hbbSettingsDict = Dict{String,Any}("verbose"=>false,"numProcesses"=>1,"conservativismLevel"=>1,
                                   "mipSettings"=>Dict{String,Any}("subsolverSettings"=>Dict{String,Any}("subsolverName"=>"OSQP")),
                                   "nlpSettings"=>Dict{String,Any}("subsolverName"=>"IPOPT"),
                                   "mipStepSettings"=>Dict{String,Any}("symbol"=>"OAmipStep","args"=>Tuple{}()),
                                   "nlpStepSettings"=>Dict{String,Any}("symbol"=>"OAnlpStep","args"=>Tuple{}()))

ssSettings = Dict{String,Any}()

# solve MIQP
workspace = OpenBB.setup("HBB",problemDict,hbbSettingsDict)

OpenBB.solveB(workspace)

solution = OpenBB.get_best_feasible_node_dict(workspace)
@assert 0.5 - workspace.settings.primalTolerance <= OpenBB.get_best_feasible_node_dict(workspace)["objUpB"] <= 0.5 + workspace.settings.primalTolerance


OpenBB.get_all_feasible_nodes_dict(workspace)


@assert OpenBB.get_numVariables(workspace) == 5
@assert OpenBB.get_numConstraints(workspace) == 2
@assert OpenBB.get_numDiscrete(workspace) == 3

@assert OpenBB.get_constraints_dependency(workspace) == [[1,2,4,5],[2,3,4]]
@assert OpenBB.get_constraint_dependency(workspace,1) == [1,2,4,5]
@assert OpenBB.get_constraint_dependency(workspace,2) == [2,3,4]
@assert OpenBB.get_objective_dependency(workspace) == [1, 2, 4, 5]

@assert OpenBB.get_variableBounds(workspace) == ([0.5, 0.0, 0.0, 0.0, -10.0], [0.5, 1.0, 1.0, 1.0, 10.0])
@assert OpenBB.get_constraintBounds(workspace) == ([0.0, 1.0], [0.0, 1.0])
@assert OpenBB.get_status_dict(workspace)["description"] == "optimalSolutionFound"


OpenBB.resetB(workspace)
@assert OpenBB.get_numFeasibleNodes(workspace) == 0
@assert OpenBB.get_status_dict(workspace)["objUpB"] == Inf

OpenBB.solveB(workspace)
OpenBB.append_constraintsB(workspace,problemDict["cnsSet"],true)
OpenBB.remove_constraintsB(workspace,[3,4],true)
OpenBB.permute_constraintsB(workspace,[2,1],true)

@assert OpenBB.get_constraintBounds(workspace) == ([1.0, 0.0], [1.0, 0.0])
OpenBB.update_boundsB(workspace,Dict("cnsLoBs"=>[0.0,0.0],"varLoBs"=>[0.5, 0.0, 0.0, 0.0, -100.0]),true)

@assert OpenBB.get_constraintBounds(workspace) == ([0.0, 0.0], [1.0, 0.0])
@assert OpenBB.get_variableBounds(workspace) == ([0.5, 0.0, 0.0, 0.0, -100.0], [0.5, 1.0, 1.0, 1.0, 10.0])
OpenBB.update_boundsB(workspace,Dict("cnsLoBs"=>[1.0,0.0]),true)
