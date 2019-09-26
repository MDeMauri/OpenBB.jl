# @Author: Massimo De Mauri <massimo>
# @Date:   2019-03-03T19:12:12+01:00
# @Email:  massimo.demauri@gmail.com
# @Filename: Gurobi_interface.jl
# @Last modified by:   massimo
# @Last modified time: 2019-09-19T12:53:46+02:00
# @License: LGPL-3.0
# @Copyright: {{copyright}}

using Gurobi
const GUROBImaxint = 2000000000

## Settings ##########################################################
mutable struct GUROBIsettings <:AbstractSettings
    # algorithm and IO
    Method::Int                    #1->simplex, 2->barrier
    OutputFlag::Bool               #Enables or disables screen output
    DisplayInterval::Int           #Frequency at which log lines are printed
    # termination
    IterationLimit::Float64        # max number of iterations allowed
    TimeLimit::Float64             # max time allowed
    # tolerances
    FeasibilityTol::Float64 	   #Primal feasibility tolerance
    OptimalityTol::Float64 	       #Dual feasibility tolerance
    PSDTol::Float64 	           #Positive semi-definite tolerance
    # simplex algorithm
    MarkowitzTol::Float64	       #Threshold pivoting tolerance
    InfUnbdInfo::Int 	           #Generate additional info for infeasible/unbounded models
    NormAdjust::Int  	           #Simplex pricing norm
    ObjScale::Float64 	           #Objective scaling
    PerturbValue::Float64 	       #Simplex perturbation magnitude
    Quad::Int 	                   #Quad precision computation in simplex
    ScaleFlag::Int 	               #Model scaling
    Sifting::Int 	               #Sifting within dual simplex
    SiftMethod::Int 	           #LP method used to solve sifting sub-problems
    SimplexPricing::Int 	       #Simplex variable pricing strategy
    # barrier algorithm
    BarIterLimit::Int 	           #Barrier iteration limit
    BarConvTol::Float64 	       #Barrier convergence tolerance
    BarQCPConvTol::Float64 	       #Barrier QCP convergence tolerance
    BarCorrectors::Int 	           #Central correction limit
    BarHomogeneous::Int 	       #Barrier homogeneous algorithm
    BarOrder::Int 	               #Barrier ordering algorithm
    Crossover::Int  	           #Barrier crossover strategy
    CrossoverBasis::Int 	       #Crossover initial basis construction strategy
    QCPDual::Int 	               #Compute dual variables for QCP models
    # thead count
    Threads::Int                   # number of threads to use
end

function GUROBIsettings(;   Method::Int=2,
                            OutputFlag::Bool=false,
                            DisplayInterval::Int=10,
                            IterationLimit::Float64=Inf,
                            TimeLimit::Float64=Inf,
                            FeasibilityTol::Float64=1e-6,
                            OptimalityTol::Float64=1e-6,
                            PSDTol::Float64=1e-6,
                            MarkowitzTol::Float64=0.0078125,
                            InfUnbdInfo::Int=0,
                            NormAdjust::Int=-1,
                            ObjScale::Float64=0.0,
                            PerturbValue::Float64=0.0002,
                            Quad::Int=-1,
                            ScaleFlag::Int=-1,
                            Sifting::Int=-1,
                            SiftMethod::Int=-1,
                            SimplexPricing::Int=-1,
                            BarIterLimit::Int=1000,
                            BarConvTol::Float64=1e-8,
                            BarQCPConvTol::Float64=1e-6,
                            BarCorrectors::Int=-1,
                            BarHomogeneous::Int=-1,
                            BarOrder::Int=-1,
                            Crossover::Int=-1,
                            CrossoverBasis::Int=0,
                            QCPDual::Int=0,
                            Threads::Int=1)::GUROBIsettings

    return GUROBIsettings(Method, OutputFlag, DisplayInterval, IterationLimit, TimeLimit,
                          FeasibilityTol, OptimalityTol, PSDTol, MarkowitzTol,
                          InfUnbdInfo, NormAdjust, ObjScale, PerturbValue, Quad,
                          ScaleFlag, Sifting, SiftMethod, SimplexPricing,
                          BarIterLimit, BarConvTol, BarQCPConvTol, BarCorrectors,
                          BarHomogeneous, BarOrder, Crossover, CrossoverBasis,
                          QCPDual, Threads)
end


## Workspace ##########################################################
# structure used for storing data for OSQP solver
mutable struct GUROBIworkspace{T1<:AbstractObjective,T2<:AbstractConstraintSet} <: AbstractWorkspace
    # problem
    problem::Problem{T1,T2}
    # memory
    environment::Gurobi.Env
    settings::GUROBIsettings
    # outdated flag
    outdated::Bool
end


## Setup & Update ##########################################################
# this function creates a Gurobi Model representing the given CvxQproblem
function setup(problem::Problem,settings::GUROBIsettings;bb_primalTolerance::Float64=Inf,bb_timeLimit=Inf)::GUROBIworkspace

    # check the problem
    @assert problem.objFun isa NullObjective || problem.objFun isa LinearObjective || problem.objFun isa QuadraticObjective
    @assert problem.cnsSet isa NullConstraintSet || problem.cnsSet isa LinearConstraintSet

    # overwrite the gurobi setting depending on the branch and bound settings
    settings.FeasibilityTol = min(settings.FeasibilityTol,bb_primalTolerance)
    if bb_timeLimit < Inf
        if settings_dict.TimeLimit == 0.
            settings_dict.TimeLimit = bb_timeLimit
        else
            settings.TimeLimit = min(settings.TimeLimit,bb_timeLimit)
        end
    end

    # reformat the settings for GUROBI
    settings_dict = Dict{Symbol,Any}()
    for field in fieldnames(GUROBIsettings)
        settings_dict[field] = getfield(settings,field)
    end

    # create the subsolver OSQPworkspace
    env = Gurobi.Env()
    Gurobi.setparams!(env;settings_dict...)

    return GUROBIworkspace(problem,env,settings,false)
end

# it marks the workspace as outdated
function make_outdated!(workspace::GUROBIworkspace)::Nothing
    workspace.outdated = true
    return
end

#
function update!(workspace::GUROBIworkspace)::Nothing
    return
end

## Solve ##########################################################
function solve!(node::BBnode,workspace::GUROBIworkspace)::Tuple{Int8,Float64}

    # update the gurobi model
    numVars = get_numVariables(workspace.problem)
    numCnss = get_numConstraints(workspace.problem)

    # create a Gurobi environment
    cnsSet = LinearConstraintSet(workspace.problem.cnsSet)
    objFun = QuadraticObjective(workspace.problem.objFun)
    model = Gurobi.gurobi_model(workspace.environment,H = objFun.Q,
                                                      f = objFun.L,
                                                      A = vcat(-cnsSet.A,cnsSet.A),
                                                      b = vcat(-node.cnsLoBs,node.cnsUpBs),
                                                      lb = node.varLoBs,
                                                      ub = node.varUpBs)

    Gurobi.update_model!(model)

    # solve problem
    runtime = @elapsed Gurobi.optimize(model)
    status = Gurobi.get_status_code(model)


    # output sol info
    if  status == 2
        status = 0 # "solved"
        node.primal = Gurobi.get_solution(model)
        node.bndDual = zeros(numVars)
        node.cnsDual = zeros(numCnss)
        node.objGap = workspace.settings.OptimalityTol
        node.objVal = (transpose(node.primal)*objFun.Q*node.primal)/2. + transpose(objFun.L)*node.primal
    elseif status in [3,4]
        status = 1 # "infeasible"
        node.primal = NaNs(numVars)
        node.bndDual = NaNs(numVars)
        node.cnsDual = NaNs(numCnss)
        node.objGap = workspace.settings.OptimalityTol
        node.objVal = Inf
    elseif status in [7,8,10,11,13]
        status = 2 # "unreliable"
        node.primal = Gurobi.get_solution(model)
        node.primal = @. min(max(node.primal,node.varLoBs),node.varUpBs)
        node.bndDual = zeros(numVars)
        node.cnsDual = zeros(numCnss)
        newObjVal = (transpose(node.primal)*objFun.Q*node.primal)/2. + transpose(objFun.L)*node.primal
        if newObjVal >= node.ObjVal - node.objGap
            node.objGap = newObjVal - node.objVal + node.objGap
            node.objVal = newObjVal
        else
            node.objGap = Inf # gurobi doesn't give enough information to estimate the gap
            node.objVal = newObjVal
            @warn "Inaccuracy in node sol, status: "*string(sol.info.status)*" (code: "*string(status)*")"
        end
    elseif status in [1,5,12]
        status = 3 # "error"
        @error "Subsover error, status: "*string(Gurobi.get_status(model))*" (code: "*string(status)*")"
    else
        @error "Subsolver unknown status: "*string(Gurobi.get_status(model))*" (code:"*string(status)*")"
    end

    return (status, runtime)
end
