# @Author: Massimo De Mauri <massimo>
# @Date:   2019-05-23T11:26:44+02:00
# @Email:  massimo.demauri@gmail.com
# @Filename: BBnode.jl
# @Last modified by:   massimo
# @Last modified time: 2021-01-13T13:44:06+01:00
# @License: LGPL-3.0
# @Copyright: {{copyright}}


# define abstract and null types
abstract type AbstractBBnode <: AbstractNode end;

# this node type is used for communication purposes
struct NullBBnode <: AbstractBBnode end

# this node type is used to arrest the branch and bound process on remote treads
struct KillerNode <:AbstractBBnode
    code::Int
end

# this node type contains the data that distinguishes a B&B problem from another
mutable struct BBnode <: AbstractBBnode
    # bounds
    varLoBs::Vector{Float}
    varUpBs::Vector{Float}
    cnsLoBs::Vector{Float}
    cnsUpBs::Vector{Float}
    # solution
    primal::Vector{Float}
    bndDual::Vector{Float}
    cnsDual::Vector{Float}
    #local cuts
    maxNumOfCuts::Int
    cutSet::LinearConstraintSet{SpMatrix{Float}}
    cutDual::Vector{Float}
    # scores
    avgFractionality::Float
    objUpB::Float
    objLoB::Float
    pseudoCost::Float
    reliable::Bool
    heuristic::Bool
    # update version
    version::Int64
end

# construct a BBnode given its lower bounds, upper bounds and solution hotstart
function BBnode(varLoBs::Vector{Float},varUpBs::Vector{Float},
                cnsLoBs::Vector{Float},cnsUpBs::Vector{Float},
                primal::Vector{Float},bndDual::Vector{Float},cnsDual::Vector{Float},
                maxNumberOfCuts::Int=0,version::Int=0)::BBnode
    # check inputs
    @assert length(varLoBs) == length(varUpBs) == length(primal) == length(bndDual)
    @assert length(cnsLoBs) == length(cnsUpBs) == length(cnsDual)
    @assert maxNumberOfCuts >= 0
    @assert version >= 0

    return BBnode(varLoBs,varUpBs,
                  cnsLoBs,cnsUpBs,
                  primal,bndDual,cnsDual,
                  maxNumberOfCuts,LinearConstraintSet(spzeros(0,length(varLoBs)),Float[],Float[]),Float[],
                  1.0,Inf,-Inf,0.0,false,false,version)
end

# overload of functions
function Base.copy(node::BBnode)::BBnode
    return BBnode(node.varLoBs,node.varUpBs,
                  node.cnsLoBs,node.cnsUpBs,
                  node.primal,node.bndDual,node.cnsDual,
                  node.maxNumOfCuts,node.cutSet,node.cutDual,
                  node.avgFractionality,node.objUpB,node.objLoB,node.pseudoCost,
                  node.reliable,node.heuristic,node.version)
end


function Base.deepcopy(node::BBnode)::BBnode
    return BBnode(copy(node.varLoBs),copy(node.varUpBs),
                  copy(node.cnsLoBs),copy(node.cnsUpBs),
                  copy(node.primal),copy(node.bndDual),copy(node.cnsDual),
                  node.maxNumOfCuts,deepcopy(node.cutSet),copy(node.cutDual),
                  node.avgFractionality,node.objUpB,node.objLoB,node.pseudoCost,
                  node.reliable,node.heuristic,node.version)
end


function add_cuts!(node::BBnode;A::SpMatrix{Float},loBs::Vector{Float},upBs::Vector{Float},forBlackList::Bool=false)::Nothing
    if get_size(node.cutSet) + size(A,1) > node.maxNumOfCuts - !forBlackList
        error("BBnode: max number of cuts violated")
    else
        node.cutSet.A = vcat(node.cutSet.A,A)
        append!(node.cutSet.loBs,loBs)
        append!(node.cutSet.upBs,upBs)
        append!(node.cutDual,0.0)
    end
    return
end





## Serialization
# this function returns the size of a serial representation of a node
function maxSerialSize_BBnode(numVars::Int,numCnss::Int,maxNumOfCuts::Int)::Int
    return 3 +                      # header
           7 +                      # average fractionality + objective + objective gap + pseudo-objective + reliable + heuristic + version
           3 +                      # num cuts + nun non-zeros in cuts + maxNumOfCuts
           4*numVars +              # variable bounds + primal + bound_dual
           3*numCnss +              # constraints bounds + constraints_dual
           3*maxNumOfCuts*numVars +   # cuts matrix (sparse)
           3*maxNumOfCuts             # cuts bounds + cuts dual

end

function maxSerialSize(node::BBnode)::Int
    return maxSerialSize_BBnode(length(node.primal),length(node.cnsDual),node.maxNumOfCuts)
end

function serialSize(node::BBnode)::Int

    # collect info
    numVars = length(node.primal)
    numCnss = length(node.cnsDual)
    numCuts = get_size(node.cutSet)
    numNZsCuts = nnz(node.cutSet.A)
    # compute size
    return 3 +                  # header
           7 +                  # average fractionality + objective + objective gap + pseudo-objective + reliable + heuristic + version
           3 +                  # num cuts + nun non-zeros in cuts + maxNumOfCuts
           4*numVars +          # variable bounds + primal + bound_dual
           3*numCnss +          # constraints bounds + constraints_dual
           3*numNZsCuts +       # cuts matrix (sparse)
           3*numCuts            # cuts bounds + cuts dual
end

function serialSize(node::NullBBnode)::Int
    return 1
end

function serialSize(node::KillerNode)::Int
    return 2
end

# this function construct and returns a SerialData object representing the node
function serialize(node::T;offset::Int=0)::SerialData where T<:AbstractBBnode

    # initialize target array
    serial = SerialData(Vector{Float}(undef,serialSize(node)+offset))
    serialize_in!(serial,node,offset=offset)
    return serial
end

# fill the destination array with the flat representation of the given node
function serialize_in!(serial::SerialData,node::BBnode;offset::Int=0)::Int

    numVars = length(node.primal)
    numCnss = length(node.cnsDual)
    numCuts = get_size(node.cutSet)
    numCutsMatrixNZs = nnz(sparse(node.cutSet.A))

    @assert length(serial) >= serialSize(node) + offset

    # header
    serial[offset+1] = 1. # type of node
    serial[offset+2] = numVars
    serial[offset+3] = numCnss
    offset += 3

    # numeric values
    serial[offset+1] = node.avgFractionality
    serial[offset+2] = node.objUpB
    serial[offset+3] = node.objLoB
    serial[offset+4] = node.pseudoCost
    serial[offset+5] = node.reliable
    serial[offset+6] = node.heuristic
    serial[offset+7] = node.version
    offset += 7

    # cuts dimensions
    serial[offset+1] = numCuts
    serial[offset+2] = nnz(sparse(node.cutSet.A))
    serial[offset+3] = node.maxNumOfCuts
    offset += 3

    # bounds
    serial[offset+1:offset+numVars] = node.varLoBs; offset+=numVars
    serial[offset+1:offset+numVars] = node.varUpBs; offset+=numVars
    serial[offset+1:offset+numCnss] = node.cnsLoBs; offset+=numCnss
    serial[offset+1:offset+numCnss] = node.cnsUpBs; offset+=numCnss

    # primal
    serial[offset+1:offset+numVars] = node.primal; offset += numVars

    # dual
    serial[offset+1:offset+numVars] = node.bndDual; offset += numVars
    serial[offset+1:offset+numCnss] = node.cnsDual; offset += numCnss

    # cuts
    cutNZs = findnz(node.cutSet.A)
    serial[offset+1:offset+numCutsMatrixNZs] = cutNZs[1]; offset += numCutsMatrixNZs
    serial[offset+1:offset+numCutsMatrixNZs] = cutNZs[2]; offset += numCutsMatrixNZs
    serial[offset+1:offset+numCutsMatrixNZs] = cutNZs[3]; offset += numCutsMatrixNZs
    serial[offset+1:offset+numCuts] = node.cutSet.loBs; offset += numCuts
    serial[offset+1:offset+numCuts] = node.cutSet.upBs; offset += numCuts
    serial[offset+1:offset+numCuts] = node.cutDual; offset += numCuts

    return offset
end

function serialize_in!(serial::SerialData,node::NullBBnode;offset::Int=0)::Int
    @assert length(serial) >= 1 + offset
    serial[offset+1] = 0. # type of node
    offset += 1
    return offset
end

function serialize_in!(serial::SerialData,node::KillerNode;offset::Int=0)::Int
    @assert length(serial) >= 2 + offset
    serial[offset+1] = -1. # type of node
    serial[offset+2] = node.code
    offset += 2
    return offset
end



# reconstruct a node from its flat representation
function BBnode(serial::SerialData;offset::Int=0)::Tuple{AbstractBBnode,Int}

    # read type
    if  serial[offset+1] == 0.0
        return (NullBBnode(),offset+1)

    elseif serial[offset+1] == 1.0

        # header
        numVars = Int(serial[offset+2])
        numCnss = Int(serial[offset+3])
        offset += 3

        # numeric values
        avgFractionality    = serial[offset+1]
        objUpB           = serial[offset+2]
        objLoB              = serial[offset+3]
        pseudoCost          = serial[offset+4]
        reliable            = Bool(serial[offset+5])
        heuristic           = Bool(serial[offset+6])
        version             = Int(serial[offset+7])
        offset += 7

        # cuts dimensions
        numCuts             = Int(serial[offset+1])
        numCutsMatrixNZs    = Int(serial[offset+2])
        maxNumOfCuts          = Int(serial[offset+3])
        offset += 3

        # bounds
        varLoBs = serial[offset+1:offset+numVars]; offset += numVars
        varUpBs = serial[offset+1:offset+numVars]; offset += numVars
        cnsLoBs = serial[offset+1:offset+numCnss]; offset += numCnss
        cnsUpBs = serial[offset+1:offset+numCnss]; offset += numCnss

        # primal
        primal = serial[offset+1:offset+numVars]; offset += numVars

        # dual
        bndDual = serial[offset+1:offset+numVars]; offset += numVars
        cnsDual = serial[offset+1:offset+numCnss]; offset += numCnss

        # cuts
        cutsNZs1 = serial[offset+1:offset+numCutsMatrixNZs]; offset += numCutsMatrixNZs
        cutsNZs2 = serial[offset+1:offset+numCutsMatrixNZs]; offset += numCutsMatrixNZs
        cutsNZs3 = serial[offset+1:offset+numCutsMatrixNZs]; offset += numCutsMatrixNZs
        cutsMatrix = sparse(cutsNZs1,cutsNZs2,cutsNZs3,numCuts,numVars)
        cutsLoBs = serial[offset+1:offset+numCuts]; offset += numCuts
        cutsUpBs = serial[offset+1:offset+numCuts]; offset += numCuts
        cutDual = serial[offset+1:offset+numCuts]; offset += numCuts

        node = BBnode(varLoBs,varUpBs,
                      cnsLoBs,cnsUpBs,
                      primal,bndDual,cnsDual,
                      maxNumOfCuts,LinearConstraintSet(cutsMatrix,cutsLoBs,cutsUpBs),cutDual,
                      avgFractionality,objUpB,objLoB,pseudoCost,
                      reliable,heuristic,version)
        return (node,offset)

    elseif serial[offset+1] == -1.0
        return (KillerNode(serial[offset+2]),offset+2)
    else
        println(serial)
        @error "type of node unknown"
        return
    end
end
