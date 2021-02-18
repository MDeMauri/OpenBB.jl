# OpenBB
An open and modular Branch and Bound framework written in Julia.
* version 1.0

## The Idea
Those are the driving ideas behind this project:
* It should be very easy to define and tackle different problem classes. Possibly, via the addition of new subsolvers.
* It must be possible to interrupt the Branch and Bound execution at any moment and restart it later.
* The state of the Branch & Bound should be transparent to the user and easy to manipulate.
* It must be possible to modify the problem (as much as it makes sense) without restarting the Branch and Bound execution.
* The code should be optimized as much as possible without harming the generality

## What is already here
* A multi-process Branch and Bound implementation for Mixed-Integer Convex/Quadratic/Linear Problems.
* A multi-process Hybrid Branch and Bound (Outer-Approximation-Like) implementation for Mixed-Integer Convex Problems.
* An integrated Python interface based on pyjulia (https://github.com/JuliaPy/pyjulia).
* OSQP binding (https://osqp.org/) based on OSQP.jl (https://github.com/oxfordcontrol/OSQP.jl).
* QPALM binding (https://github.com/Benny44/QPALM) based on QPALM.jl (https://github.com/kul-forbes/QPALM.jl).
* CLP binding (https://projects.coin-or.org/Clp) modified from CLP.jl(https://github.com/xhub/Clp.jl)
* Ipopt binding (https://github.com/coin-or/Ipopt) based on Ipopt.jl (https://github.com/jump-dev/Ipopt.jl)
* CLP (https://projects.coin-or.org/Clp) and CPLEX (https://www.ibm.com/analytics/cplex-optimizer) interfaces via JuMP (https://github.com/jump-dev/JuMP.jl) based on Clp.jl (https://github.com/xhub/Clp.jl) and CPLEX.jl (https://github.com/jump-dev/CPLEX.jl)


The code allows you to easily define your custom priority rules for the selection of the next node to solve and the next variable to branch on. Moreover, it is very easy to define custom stopping criteria for Branch & Bound. Further, the whole code is subsolver-invariant and a new subsolver can be added by simply overloading the interface functions that the user can find in the "subsolvers_interfaces" folder. Finally, OpenBB provides a number of functions that allow safe constraints addition, problem expansion, bounds restrictions, etc..



### Disclaimer 1
The project is a early phase of development. Some features are still missing. Here is a tentative to do list:
* Add Semi-Definite and Second-Order-Cone Programming support
* Add bound propagation, preprocessing and cuts generation techniques
* Document all

Feel free to suggest missing functionalities and to collaborate in the project.

### Disclaimer 2
If you are going to use OpenBB for your reasearch please cite me: Massimo De Mauri (massimo.demauri@protonmail.be). (and of course, please, cite the authors of the subsolvers you are going to use)
