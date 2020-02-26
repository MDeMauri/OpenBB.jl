# @Author: Massimo De Mauri <massimo>
# @Date:   2019-05-23T11:26:44+02:00
# @Email:  massimo.demauri@gmail.com
# @Filename: BBnode.jl
# @Last modified by:   massimo
# @Last modified time: 2020-02-26T21:42:00+01:00
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




# Serialization
# this function returns the size of a serial representation of a node
function max_serial_size_BBnode(numVars::Int,numCnss::Int,numCuts::Int)::Int
    return 3 +                  # header
           6 +                  # average fractionality + objective + objective gap + pseudo-objective + reliable + version
           2 +                  # num cuts + nun non-zeros in cuts
           4*numVars +          # variable bounds + primal + bound_dual
           3*numCnss +          # constraints bounds + constraints_dual
           3*numCuts*numVars +  # cuts matrix (sparse)
           3*numCuts            # cuts bounds + cuts dual

end

function max_serial_size(node::BBnode)::Int
    return max_serial_size_BBnode(length(node.primal),length(node.cnsDual),get_size(node.cuts))
end

function serial_size(node::BBnode)::Int

    # collect info
    numVars = length(node.primal)
    numCnss = length(node.cnsDual)
    numCuts = get_size(node.cuts)
    numNZsCuts = nnz(node.cuts.A)
    # compute size
    return 3 +                  # header
           6 +                  # average fractionality + objective + objective gap + pseudo-objective + reliable + version
           2 +                  # num cuts + nun non-zeros in cuts
           4*numVars +          # variable bounds + primal + bound_dual
           3*numCnss +          # constraints bounds + constraints_dual
           3*numNZsCuts +       # cuts matrix (sparse)
           3*numCuts            # cuts bounds + cuts dual
end

function serial_size(node::NullBBnode)::Int
    return 1
end

function serial_size(node::KillerNode)::Int
    return 2
end

# this function construct and returns a SerialData object representing the node
function serialize(node::T;offset::Int=0)::SerialData where T<:AbstractBBnode

    # initialize target array
    serial = SerialData(Array{Float64,1}(undef,serial_size(node)+offset))
    serialize_in!(serial,node,offset=offset)
    return serial
end

# fill the destination array with the flat representation of the given node
function serialize_in!(serial::SerialData,node::BBnode;offset::Int=0)::Int

    numVars = length(node.primal)
    numCnss = length(node.cnsDual)
    numCuts = get_size(node.cuts)
    numCutsMatrixNZs = nnz(sparse(node.cuts.A))

    @assert length(serial) >= serial_size(node) + offset

    # header
    serial[offset+1] = 1. # type of node
    serial[offset+2] = numVars
    serial[offset+3] = numCnss
    offset += 3

    # numeric values
    serial[offset+1] = node.avgAbsFrac
    serial[offset+2] = node.objVal
    serial[offset+3] = node.objGap
    serial[offset+4] = node.pseudoObjective
    serial[offset+5] = node.reliable
    serial[offset+6] = node.version
    offset += 6

    # cuts dimensions
    serial[offset+1] = numCuts
    serial[offset+2] = nnz(sparse(node.cuts.A))
    offset += 2

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
    cutNZs = findnz(node.cuts.A)
    serial[offset+1:offset+numCutsMatrixNZs] = cutNZs[1]; offset += numCutsMatrixNZs
    serial[offset+1:offset+numCutsMatrixNZs] = cutNZs[2]; offset += numCutsMatrixNZs
    serial[offset+1:offset+numCutsMatrixNZs] = cutNZs[3]; offset += numCutsMatrixNZs
    serial[offset+1:offset+numCuts] = node.cuts.loBs; offset += numCuts
    serial[offset+1:offset+numCuts] = node.cuts.upBs; offset += numCuts
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
    serial[offset+2] = node.count
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
        avgAbsFrac          = serial[offset+1]
        objective           = serial[offset+2]
        objGap              = serial[offset+3]
        pseudoObjective     = serial[offset+4]
        reliable            = Bool(serial[offset+5])
        version             = Int(serial[offset+6])
        offset += 6

        # cuts dimensions
        numCuts             = Int(serial[offset+1])
        numCutsMatrixNZs    = Int(serial[offset+2])
        offset += 2

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
                      LinearConstraintSet(cutsMatrix,cutsLoBs,cutsUpBs),cutDual,
                      avgAbsFrac,objective,objGap,pseudoObjective,
                      reliable,version)
        return (node,offset)

    elseif serial[offset+1] == -1.0
        return (KillerNode(serial[offset+2]),offset+2)
    else
        println(serial)
        @error "type of node unknown"
        return
    end
end
