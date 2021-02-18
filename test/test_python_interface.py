# @Author: Massimo De Mauri <massimo>
# @Date:   2019-06-19T18:15:50+02:00
# @Email:  massimo.demauri@gmail.com
# @Filename: test_interface.py
# @Last modified by:   massimo
# @Last modified time: 2020-12-29T11:07:59+01:00
# @License: LGPL-3.0

import numpy as np
from python_interface import *


# load openBB
OpenBB = OpenBBinterface()


# create the problem
problem = {"objFun":{},"cnsSet":{},"varSet":{}}

problem["objFun"]["Q"] = 2*np.eye(5)
problem["objFun"]["Q"][2,2] = 0
problem["objFun"]["L"] = np.zeros(5)

problem["cnsSet"]["A"] = np.matrix([[1.,1.,0.,-1.,-1.],[0.,1.,1.,1.,0.]])
problem["cnsSet"]["loBs"] = np.array([0.,1.])
problem["cnsSet"]["upBs"] = np.array([0.,1.])

problem["varSet"]["loBs"] = [.5,0,0,0,-10]
problem["varSet"]["upBs"] = [.5,1,1,1,10]
problem["varSet"]["vals"] = [.5,0,1,0,.5]
problem["varSet"]["dscIndices"] = [1,2,3]
problem["varSet"]["sos1Groups"] = [1,1,1]

settings = {"subsolver":"osqp",
            "verbose":False,
            "numProcesses":1,
            "conservativismLevel":1
            }

OpenBB.setup(problem,settings)
OpenBB.solve()

assert OpenBB.get_best_feasible_node()["objective"] >= 0.5 - OpenBB.get_settings()["primalTolerance"]
assert OpenBB.get_best_feasible_node()["objective"] <= 0.5 + OpenBB.get_settings()["primalTolerance"]

OpenBB.get_all_feasible_nodes()
OpenBB.get_best_node()

assert OpenBB.get_numVariables() == 5
assert OpenBB.get_numConstraints() == 2
assert OpenBB.get_numDiscrete() == 3

OpenBB.get_constraints_sparsity()
assert all(OpenBB.get_constraint_dependency(0) == [0,1,3,4])
assert all(OpenBB.get_constraint_dependency(1) == [1,2,3])
assert all(OpenBB.get_objective_dependency()[0] == [0, 1, 3, 4])

assert all(OpenBB.get_variableBounds()[0] == [0.5, 0.0, 0.0, 0.0, -10.0])
assert all(OpenBB.get_constraintBounds()[0] == [0.0, 1.0])
assert OpenBB.get_status()["description"] == "optimalSolutionFound"

OpenBB.reset_explored_nodes()
assert OpenBB.get_numInactiveNodes() == OpenBB.get_numSolutions() == 0
assert OpenBB.get_status()["objUpB"] == np.inf

OpenBB.clear()
assert OpenBB.get_numActiveNodes() == 0

OpenBB.reset()
assert OpenBB.get_numActiveNodes() == 1

OpenBB.solve()
OpenBB.append_constraints(problem["cnsSet"],False,True,False)
OpenBB.remove_constraints([2,3],True,True,False)
OpenBB.permute_constraints([1,0],True,True,False)
assert all(OpenBB.get_constraintBounds()[0] == [1.0, 0.0])


OpenBB.update_bounds({"cnsLoBs":[0.0,0.0],"varLoBs":[0.5, 0.0, 0.0, 0.0, -100.0]},False,True,False)
assert all(OpenBB.get_constraintBounds()[0] == [0.0, 0.0])
assert all(OpenBB.get_variableBounds()[0] == [0.5, 0.0, 0.0, 0.0, -100.0])
OpenBB.update_bounds({"cnsLoBs":[1.0,0.0]},False,True,False)


OpenBB.append_problem(problem,False,False,False)
OpenBB.solve()
assert OpenBB.get_best_feasible_node()["objective"] >= 1.0 - OpenBB.get_settings()["primalTolerance"]
assert OpenBB.get_best_feasible_node()["objective"] <= 1.0 + OpenBB.get_settings()["primalTolerance"]


OpenBB.integralize_variables([4])
OpenBB.solve()
assert OpenBB.get_status()["description"] == "infeasible"


OpenBB.get_constraints()
OpenBB.get_objective()

print(" - Python interface, ok")
