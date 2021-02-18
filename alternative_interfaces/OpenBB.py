# @Author: Massimo De Mauri <massimo>
# @Date:   2019-06-19T16:35:26+02:00
# @Email:  massimo.demauri@gmail.com
# @Filename: python_interface.py
# @Last modified by:   massimo
# @Last modified time: 2021-02-01T13:01:52+01:00
# @License: LGPL-3.0
# @Copyright: {{copyright}}


# import the needed components
from numpy import array, matrix
from julia import OpenBB

class OpenBBinterface:

    def __init__(self):
        self.jl = OpenBB
        self.jl.include_flat_interface()
        self.workspace = None


    ######################## general ########################

    def eval_string(self,string):
        return self.jl.eval_string(string)


    ######################## setup ########################
    # load subsolver interface
    def load_subsolver_interface(self,subsolver_name):
        self.jl.load_subsolver_interface(subsolver_name)
        return


    # setup a OpenBB workspace
    def setup(self,algorithm,problemDict,options={}):

        # clear memory
        self.jl.unload_all_libs()

        # load the necessary subsolver interfaces
        if algorithm == "bb" or algorithm == "BB":
            self.load_subsolver_interface(options["subsolverSettings"]["subsolverName"])
        elif algorithm == "hbb" or algorithm == "HBB":
            self.load_subsolver_interface(options["nlpSettings"]["subsolverName"])
            self.load_subsolver_interface(options["mipSettings"]["subsolverSettings"]["subsolverName"])
        else:
            raise NameError("OpenBB: algorithm unknown.")

        # setup new workspace
        self.workspace = self.jl.setup(algorithm,problemDict,options)

        return


    ######################## solve ########################
    # solve the problem
    def solve(self,initial_point=[]):
        if self.workspace is None:
            raise NameError("workspace not initialized, please run setup(algorithm,problemDict,options)")

        # solve the problem
        self.jl.solveB(self.workspace,array(initial_point).flatten())
        return



    ######################## inspect ########################

    # get statistics and results of the last optimization
    def get_status(self):
        return self.jl.get_status_dict(self.workspace)

    # ...
    def get_best_feasible_node(self):
        return self.jl.get_best_feasible_node_dict(self.workspace)

    # ...
    def get_all_feasible_nodes(self):
        return self.jl.get_all_feasible_nodes_dict(self.workspace)

    # ...
    def get_best_node(self):
        return self.jl.get_best_node_dict(self.workspace)

    # ...
    def get_numVariables(self):
        return self.jl.get_numVariables(self.workspace)

    # ...
    def get_numConstraints(self):
        return self.jl.get_numConstraints(self.workspace)

    # ...
    def get_numDiscrete(self):
        return self.jl.get_numDiscrete(self.workspace)


    # ...
    def get_constraints(self):
        return self.jl.get_constraints_dict(self.workspace)

    # ...
    def get_objective(self):
        return self.jl.get_objective_dict(self.workspace)


    # ...
    def get_constraint_dependency(self,index):
        out = self.jl.get_constraint_dependency(self.workspace,index+1)
        for k in range(len(out)):
            out[k] = out[k] - 1
        return out

    # ...
    def get_objective_dependency(self):
        out = self.jl.get_objective_dependency(self.workspace)
        for k in range(len(out[0])):
            out[0][k] = out[0][k] - 1
            out[1][k] = out[1][k] - 1
        return out

    # ...
    def get_variableBounds(self):
        return self.jl.get_variableBounds(self.workspace)

    # ...
    def get_constraintBounds(self):
        return self.jl.get_constraintBounds(self.workspace)

    # ...
    def get_numActiveNodes(self):
        return self.jl.get_numActiveNodes(self.workspace)

    # ...
    def get_numInactiveNodes(self):
        return self.jl.get_numInactiveNodes(self.workspace)

    # ...
    def get_numSolutions(self):
        return self.jl.get_numSolutions(self.workspace)


    ######################## update workspace ########################

    # ...
    def update(self):
        self.jl.updateB(self.workspace)
        return

    # ...
    def reset(self):
        self.jl.resetB(self.workspace)
        return

    # ...
    def clear(self):
        self.jl.clearB(self.workspace)
        return


    ######################## update problem ########################
    # ...
    def append_constraints(self,newConstraintsDict):
        if not 'A' in newConstraintsDict or not 'loBs' in newConstraintsDict or not 'upBs' in newConstraintsDict:
            raise NameError('newConstraintsDict has to be a dictionary with the following keywords: A,loBs,upBs')

        # reformat the constraints set
        if len(newConstraintsDict)>0:
            localConstraintsDict = {}
            localConstraintsDict['A'] = matrix(newConstraintsDict['A'])
            localConstraintsDict['loBs'] = array(newConstraintsDict['loBs']).flatten()
            localConstraintsDict['upBs'] = array(newConstraintsDict['upBs']).flatten()

        self.jl.append_constraintsB(self.workspace,localConstraintsDict)
        return

    # ...
    def insert_constraints(self,newConstraintsDict,index):
        # check correctness of the input
        if not 'A' in newConstraintsDict or not 'loBs' in newConstraintsDict or not 'upBs' in newConstraintsDict:
            raise NameError('newConstraintsDict has to be a dictionary with the following keywords: A,loBs,upBs')

        # reformat the constraints set
        if len(newConstraintsDict)>0:
            localConstraintsDict = {}
            localConstraintsDict['A'] = matrix(newConstraintsDict['A'])
            localConstraintsDict['loBs'] = array(newConstraintsDict['loBs']).flatten()
            localConstraintsDict['upBs'] = array(newConstraintsDict['upBs']).flatten()

        self.jl.insert_constraintsB(self.workspace,localConstraintsDict,index+1)
        return


    # ...
    def remove_constraints(self,indices):
        for k in range(len(indices)):
            indices[k] = indices[k] + 1
        self.jl.remove_constraintsB(self.workspace,indices)
        return

    # ...
    def permute_constraints(self,permutation):
        for k in range(len(permutation)):
            permutation[k] = permutation[k] + 1
        self.jl.permute_constraintsB(self.workspace,permutation)
        return

    # ...
    def update_bounds(self,boundsDict):
        self.jl.update_boundsB(self.workspace,boundsDict)
        return

    # ...
    def set_objective(self,newObjectiveDict):
        # reformat the objective info
        if len(newObjectiveDict) > 0:
            localObjectiveDict = {}
            if 'Q' in newObjectiveDict: localObjectiveDict['Q'] = matrix(newObjectiveDict['Q'])
            if 'L' in newObjectiveDict: localObjectiveDict['L'] = array(newObjectiveDict['L']).flatten()

            self.jl.set_objectiveB(self.workspace,localObjectiveDict)
        return

    # ...
    def set_constraintSet(self,newConstraintsDict):
        # reformat the constraints set
        if len(newConstraintsDict)>0:
            localConstraintsDict = {}
            localConstraintsDict['A'] = matrix(newConstraintsDict['A'])
            localConstraintsDict['loBs'] = array(newConstraintsDict['loBs']).flatten()
            localConstraintsDict['upBs'] = array(newConstraintsDict['upBs']).flatten()

            self.jl.set_constraintSetB(self.workspace,localConstraintsDict)
        return


    # adds new variables, new constraints and a new part of the objective to the problem
    def append_problem(self,problemDict):
        if not problemDict is None:
            # check the input
            if not 'objFun' in problemDict:
                raise NameError('keyword \'objFun\' not found')
            elif not 'cnsSet' in problemDict:
                raise NameError('keyword \'cnsSet\' not found')
            elif not 'varSet' in problemDict:
                raise NameError('keyword \'varSet\' not found')
            else:

                # reformat the objective
                localProblemDict = {'objFun':{},'cnsSet':{},'varSet':{}}
                if 'Q' in problemDict['objFun']: localProblemDict['objFun']['Q'] = matrix(problemDict['objFun']['Q'])
                localProblemDict['objFun']['L'] = array(problemDict['objFun']['L']).flatten()

                # reformat the constraints set
                localProblemDict['cnsSet']['A'] = matrix(problemDict['cnsSet']['A'])
                localProblemDict['cnsSet']['loBs'] = array(problemDict['cnsSet']['loBs']).flatten()
                localProblemDict['cnsSet']['upBs'] = array(problemDict['cnsSet']['upBs']).flatten()

                # reformat the variables set
                localProblemDict['varSet']['loBs'] = array(problemDict['varSet']['loBs']).flatten()
                localProblemDict['varSet']['upBs'] = array(problemDict['varSet']['upBs']).flatten()
                if 'vals' in problemDict['varSet']: localProblemDict['varSet']['vals'] = array(problemDict['varSet']['vals']).flatten()
                if 'dscIndices' in problemDict['varSet']: localProblemDict['varSet']['dscIndices'] = [int(problemDict['varSet']['dscIndices'][k])+1 for k in range(len(problemDict['varSet']['dscIndices']))]
                if 'sos1Groups' in problemDict['varSet']: localProblemDict['varSet']['sos1Groups'] = [int(problemDict['varSet']['sos1Groups'][k]) for k in range(len(problemDict['varSet']['sos1Groups']))]
                if 'pseudoCosts' in problemDict['varSet']: localProblemDict['varSet']['pseudoCosts'] = array(problemDict['varSet']['pseudoCosts']).flatten()

                self.jl.append_problemB(self.workspace,localProblemDict)
        return


    # marks as discrete formally not discrete variables
    def integralize_variables(self,newDscIndices,newSos1Groups=array([],int)):
        for k in range(len(newDscIndices)):
            newDscIndices[k] = newDscIndices[k] + 1
        self.jl.integralize_variablesB(self.workspace,newDscIndices,newSos1Groups)
        return

    # ...
    def update_objectiveCutoff(self,newCutoff):
        self.jl.update_objectiveCutoffB(self.workspace,newCutoff)
        return
