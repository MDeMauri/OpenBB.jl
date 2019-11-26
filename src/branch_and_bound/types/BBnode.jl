# @Author: Massimo De Mauri <massimo>
# @Date:   2019-05-23T11:26:44+02:00
# @Email:  massimo.demauri@gmail.com
# @Filename: BBnode.jl
# @Last modified by:   massimo
# @Last modified time: 2019-11-22T14:43:11+01:00
# @License: LGPL-3.0
# @Copyright: {{copyright}}


# define abstract and null types
abstract type AbstractBBnode end; struct NullBBnode <: AbstractBBnode end

# this is the set of data that distinguishes a node from another
mutable struct BBnode <: AbstractBBnode
    # bounds
    varLoBs::Array{Float64,1}
    varUpBs::Array{Float64,1}
    cnsLoBs::Array{Float64,1}
    cnsUpBs::Array{Float64,1}
    # solution
    primal::Array{Float64,1}
    bndDual::Array{Float64,1}
    cnsDual::Array{Float64,1}
    #local cuts
    cuts::LinearConstraintSet{SparseMatrixCSC{Float64,Int}}
    cutDual::Array{Float64,1}
    # scores
    avgAbsFrac::Float64
    objVal::Float64
    objGap::Float64
    pseudoObjective::Float64
    reliable::Bool
    # update version
    version::Int64
end

# construct a BBnode given its lower bounds, upper bounds and solution hotstart
function BBnode(varLoBs::Array{Float64,1},varUpBs::Array{Float64,1},
                cnsLoBs::Array{Float64,1},cnsUpBs::Array{Float64,1},
                primal::Array{Float64,1},bndDual::Array{Float64,1},cnsDual::Array{Float64,1},
                maxNumberOfCuts::Int,version::Int64)::BBnode
    # check inputs
    @assert length(varLoBs) == length(varUpBs) == length(primal) == length(bndDual)
    @assert length(cnsLoBs) == length(cnsUpBs) == length(cnsDual)
    @assert maxNumberOfCuts >= 0
    @assert version >= 0

    return BBnode(varLoBs,varUpBs,
                  cnsLoBs,cnsUpBs,
                  primal,bndDual,cnsDual,
                  LinearConstraintSet(spzeros(maxNumberOfCuts,length(varLoBs)),zeros(maxNumberOfCuts),zeros(maxNumberOfCuts)),zeros(maxNumberOfCuts),
                  NaN,NaN,0.0,NaN,
                  false,version)
end

# this node is used to arrest the branch and bound process
struct KillerNode <:AbstractBBnode
    count::Int
end

# overload of functions
import Base.copy
function copy(node::BBnode)::BBnode
    return BBnode(node.varLoBs,node.varUpBs,
                  node.cnsLoBs,node.cnsUpBs,
                  node.primal,node.bndDual,node.cnsDual,
                  node.cuts,node.cutDual,
                  node.avgAbsFrac,node.objVal,node.pseudoObjective,
                  node.reliable,node.version)
end


import Base.deepcopy
function deepcopy(node::BBnode)::BBnode
    return BBnode(copy(node.varLoBs),copy(node.varUpBs),
                  copy(node.cnsLoBs),copy(node.cnsUpBs),
                  copy(node.primal),copy(node.bndDual),copy(node.cnsDual),
                  deepcopy(node.cuts),copy(node.cutDual),
                  node.avgAbsFrac,node.objVal,node.objGap,node.pseudoObjective,
                  node.reliable,node.version)
end
