# @Author: Massimo De Mauri <massimo>
# @Date:   2019-06-19T16:35:26+02:00
# @Email:  massimo.demauri@gmail.com
# @Filename: python_interface.py
# @Last modified by:   massimo
# @Last modified time: 2021-01-27T00:47:34+01:00
# @License: LGPL-3.0
# @Copyright: {{copyright}}


# import the needed components
from os import path
from copy import copy
from warnings import warn
from numpy import array, matrix
from julia import Julia,OpenBB


class MPC_addon:

    def __init__(self,openbb_interface):
        self.openbb = openbb_interface
        return

    def BB_mpc_shift(self,shift_steps,tail_costs_dict,tail_constraints_dict,
                     reference_solution=array([],float),measured_state=array([],float),
                     mode="fullRTI",subsolver_iterations=1,suppress_warnings=False):

        # reformat the reference solution and the measured state
        reference_solution_ = array(reference_solution).flatten()
        measured_state_ = array(measured_state).flatten()

        self.openbb.jl.mpc_shiftB(self.openbb.workspace,shift_steps,tail_costs_dict,tail_constraints_dict,
                                  reference_solution_,measured_state_,
                                  mode,subsolver_iterations,suppress_warnings)
        return



    def HBB_mpc_shift_assisted(self,shift_steps,shifted_problem_dict,
                               reference_solution=array([],float),measured_state=array([],float),
                               mode="fullRTI",subsolver_iterations=1,
                               suppress_warnings=False):

        # reformat the reference solution and the measured state
        reference_solution_ = array(reference_solution).flatten()
        measured_state_ = array(measured_state).flatten()

        self.openbb.jl.mpc_shiftB(self.openbb.workspace,shift_steps,shifted_problem_dict,
                                  reference_solution_,measured_state_,
                                  mode,subsolver_iterations,
                                  suppress_warnings)
        return


    def oc_insert_constraints(self,newConstraintsDict,numVarsPerStep,
                              suppress_warnings=False,local_only=False):

        # reformat the constraints set
        local_constraints_dict = {}
        local_constraints_dict['A'] = matrix(newConstraintsDict['A'])
        local_constraints_dict['loBs'] = array(newConstraintsDict['loBs']).flatten()
        local_constraints_dict['upBs'] = array(newConstraintsDict['upBs']).flatten()

        self.openbb.jl.insert_constraints_oc(local_constraints_dict,numVarsPerStep,suppress_warnings,local_only)
        return
