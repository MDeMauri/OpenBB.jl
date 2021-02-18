# @Author: Massimo De Mauri <massimo>
# @Date:   2020-02-26T21:24:43+01:00
# @Email:  massimo.demauri@gmail.com
# @Filename: test_serialization.jl
# @Last modified by:   massimo
# @Last modified time: 2020-12-28T17:32:13+01:00
# @License: LGPL-3.0
# @Copyright: {{copyright}}

using OpenBB
using SparseArrays

# test problem serialization
varSet = [OpenBB.VariableSet(loBs = [-Inf, 0.0, 0.0], upBs = [Inf, 5.0, 1.0],dscIndices = [2, 3])]

cnsMat = [0.5325601466335066 0.9869379665568074 0.42237296990956597; 0.47711227338871587 0.2068112421256656 0.9385725471537107]
cnsLoBs = [-0.13427182108222568, -0.4968347607289445]
cnsUpBs = [1.2664268833269365, 1.463593794258908]

cnsSet = [OpenBB.LinearConstraintSet(A=cnsMat,loBs=cnsLoBs,upBs=cnsUpBs),
          OpenBB.LinearConstraintSet(A=sparse(cnsMat),loBs=cnsLoBs,upBs=cnsUpBs)]


qTerm = [0.5028313177754078 0.0 0.0; 0.0 0.8881848812969096 0.0; 0.0 0.0 0.3593351390040926]
lTerm = [0.7324805425775212, 0.7804441548831766, 0.46369934900873333]
objFun = [OpenBB.LinearObjective(L=lTerm),OpenBB.LinearObjective(L=sparsevec(lTerm)),
          OpenBB.QuadraticObjective(Q=qTerm,L=lTerm),OpenBB.QuadraticObjective(Q=sparse(qTerm),L=sparsevec(lTerm))]

for var in varSet
    for cns in cnsSet
        for obj in objFun
            orig = OpenBB.Problem(varSet=var,cnsSet=cns,objFun=obj)
            res,_ = OpenBB.Problem(OpenBB.serialize(orig))



            for field in fieldnames(typeof(var))
               if getfield(orig.varSet,field) != getfield(res.varSet,field)
                   @error "serialization error: "*string(typeof(var))*" mismatch in "*string(field)
               end
            end

            for field in fieldnames(typeof(cns))
               if getfield(orig.cnsSet,field) != getfield(res.cnsSet,field)
                   @error "serialization error: "*string(typeof(cns))*" mismatch in "*string(field)
               end
            end

            for field in fieldnames(typeof(obj))
               if getfield(orig.objFun,field) != getfield(res.objFun,field)
                   @error "serialization error: "*string(typeof(obj))*" mismatch in "*string(field)
               end
            end
        end
    end
end


# nodes serialization
node1 = OpenBB.NullBBnode()
node2,_ = OpenBB.BBnode(OpenBB.serialize(node1))
for field in fieldnames(OpenBB.NullBBnode)
    if getfield(node1,field) != getfield(node2,field)
        @error "serialization error: OpenBB.NullBBnode mismatch in "*string(field)
    end
end

node1 = OpenBB.KillerNode(2)
node2,_ = OpenBB.BBnode(OpenBB.serialize(node1))
for field in fieldnames(OpenBB.KillerNode)
    if getfield(node1,field) != getfield(node2,field)
        @error "serialization error: OpenBB.KillerNode mismatch in "*string(field)
    end
end

node1 = OpenBB.BBnode([1.,2.,3.],[4.,5.,6.],[7.,8.],[9.,10.],
                      [11.,12.,13.],[14.,15.,16.],[17.,18.],
                      5,cnsSet[2],[19.,20.],.21,22.,23.,24.,true,false,0)
node2,_ = OpenBB.BBnode(OpenBB.serialize(node1))
for field in fieldnames(OpenBB.BBnode)
    if field == :cutSet
        for field in fieldnames(OpenBB.LinearConstraintSet)
            if getfield(node1.cutSet,field) != getfield(node2.cutSet,field)
                @error "serialization error: OpenBB.BBnode mismatch in cutSet."*string(field)
            end
        end
    else
        if getfield(node1,field) != getfield(node2,field)
            @info getfield(node1,field)
            @info getfield(node2,field)
            @error "serialization error: OpenBB.BBnode mismatch in "*string(field)
        end
    end
end
