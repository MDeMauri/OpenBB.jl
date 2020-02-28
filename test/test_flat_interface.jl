# @Author: Massimo De Mauri <massimo>
# @Date:   2019-06-19T18:15:50+02:00
# @Email:  massimo.demauri@gmail.com
# @Filename: test_flat_interface.jl
# @Last modified by:   massimo
# @Last modified time: 2019-09-25T17:37:00+02:00
# @License: LGPL-3.0

using OpenBB
using LinearAlgebra
using SparseArrays
using OpenBB


problem = Dict("objFun"=>Dict(),"cnsSet"=>Dict(),"varSet"=>Dict())

problem["objFun"]["Q"] = 2. *OpenBB.speye(5)
problem["objFun"]["Q"][3,3] = 0
dropzeros!(problem["objFun"]["Q"])
problem["objFun"]["L"] = zeros(5)

problem["cnsSet"]["A"] = vcat([1. 1. 0. -1. -1.],[0. 1. 1. 1. 0.])
problem["cnsSet"]["loBs"] = [0.,1.]
problem["cnsSet"]["upBs"] = [0.,1.]

problem["varSet"]["loBs"] = [.5,0,0,0,-10]
problem["varSet"]["upBs"] = [.5,1,1,1,10]
problem["varSet"]["vals"] = [.5,0,1,0,.5]
problem["varSet"]["dscIndices"] = [2,3,4]
problem["varSet"]["sos1Groups"] = [1,1,1]

bbSettings = Dict("verbose"=>false,"numProcesses"=>1,"interactiveMode"=>true)
ssSettings = Dict()


if OpenBB.withOSQP()
    subsolver = "osqp"
elseif OpenBB.withQPALM()
    subsolver = "qpalm"
elseif OpenBB.withGUROBI()
    subsolver = "gurobi"
end

OpenBB.setup("osqp",problem,bbSettings,ssSettings)
OpenBB.solve_b()


@assert 0.5 - OpenBB.get_settings()["primalTolerance"] <= OpenBB.get_best_solution()["objUpB"] <= 0.5 + OpenBB.get_settings()["primalTolerance"]
OpenBB.get_subsolver_settings()

OpenBB.get_all_solutions()
OpenBB.get_best_node()

@assert OpenBB.get_numVariables() == 5
@assert OpenBB.get_numConstraints() == 2
@assert OpenBB.get_numDiscreteVariables() == 3

OpenBB.get_constraints_sparsity()
@assert OpenBB.get_constraint_sparsity(1) == [1,2,4,5]
@assert OpenBB.get_constraint_sparsity(2) == [2,3,4]
@assert OpenBB.get_objective_sparsity() == (([1, 2, 4, 5], [1, 2, 4, 5]),Int64[])

@assert OpenBB.get_variableBounds() == ([0.5, 0.0, 0.0, 0.0, -10.0], [0.5, 1.0, 1.0, 1.0, 10.0])
@assert OpenBB.get_constraintBounds() == ([0.0, 1.0], [0.0, 1.0])
@assert OpenBB.get_status()["description"] == "optimalSolutionFound"

OpenBB.reset_b()
@assert OpenBB.get_numUnactiveNodes() == OpenBB.get_numSolutions() == 0
@assert OpenBB.get_status()["objUpB"] == Inf

OpenBB.clear_b()
@assert OpenBB.get_numActiveNodes() == 0

OpenBB.reset_b()
@assert OpenBB.get_numActiveNodes() == 1

OpenBB.solve_b()
OpenBB.append_constraints_b(problem["cnsSet"],true,false)
OpenBB.remove_constraints_b([3,4],true,false)
OpenBB.permute_constraints_b([2,1],true,false)
@assert OpenBB.get_constraintBounds() == ([1.0, 0.0], [1.0, 0.0])
OpenBB.update_bounds_b(Dict("cnsLoBs"=>[0.0,0.0],"varLoBs"=>[0.5, 0.0, 0.0, 0.0, -100.0]),true,false)
@assert OpenBB.get_constraintBounds() == ([0.0, 0.0], [1.0, 0.0])
@assert OpenBB.get_variableBounds() == ([0.5, 0.0, 0.0, 0.0, -100.0], [0.5, 1.0, 1.0, 1.0, 10.0])
OpenBB.update_bounds_b(Dict("cnsLoBs"=>[1.0,0.0]),true,false)


OpenBB.append_problem_b(problem,true,false)
OpenBB.update_b()
OpenBB.solve_b()
@assert 1.0 - OpenBB.get_settings()["primalTolerance"] <= OpenBB.get_best_solution()["objUpB"] <= 1.0 + OpenBB.get_settings()["primalTolerance"]

OpenBB.integralize_variables_b([5])
OpenBB.solve_b()
@assert OpenBB.get_status()["description"] == "infeasible"

OpenBB.get_constraints()
OpenBB.get_objective()

println(" - flat interface, ok")
