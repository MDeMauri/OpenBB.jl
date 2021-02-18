# @Author: Massimo De Mauri <massimo>
# @Date:   2019-03-05T15:04:00+01:00
# @Email:  massimo.demauri@gmail.com
# @Filename: test_problem_definitions_fundamentals.jl
# @Last modified by:   massimo
# @Last modified time: 2020-10-30T22:51:47+01:00
# @License: LGPL-3.0
# @Copyright: {{copyright}}

using OpenBB
using LinearAlgebra
using SparseArrays

obj = Array{OpenBB.ObjectiveFunction,1}(undef,2)
obj[1] = OpenBB.QuadraticObjective(Q=Matrix(1.0I,4,4,),L=ones(4))
obj[2] = OpenBB.LinearObjective(L=ones(4))

cns = Array{OpenBB.ConstraintSet,1}(undef,1)
cns[1] = OpenBB.LinearConstraintSet(A=sparse(hcat(ones(1,2),zeros(1,2))),loBs=[1.],upBs=[1.])

var = Array{OpenBB.VariableSet,1}(undef,1)
var[1] = OpenBB.VariableSet(loBs=zeros(4),upBs=Infs(4),vals=zeros(4),dscIndices=[1,3],sos1Groups=[0,1])


# test fundamentals
for v in var
    tmp = copy(deepcopy(v))
    typeof(tmp)(tmp)
    @assert OpenBB.get_size(tmp) == 4
    @assert OpenBB.get_numDiscrete(tmp) == 2
    @assert OpenBB.get_bounds(tmp) == (zeros(4),Infs(4))
    @assert OpenBB.get_discreteIndices(tmp) == [1,3]
    @assert OpenBB.get_sos1Groups(tmp) == [0,1]
    @assert OpenBB.get_pseudoCosts(tmp) == (1e-4*ones(2,2),Int.(zeros(2,2)))
end

for c in cns
    tmp = copy(deepcopy(c))
    tmp = sparse(tmp)
    @assert OpenBB.get_size(tmp) == 1
    @assert OpenBB.get_numVariables(tmp) == 4
    @assert OpenBB.get_bounds(tmp) ==  ([1.],[1.])
    if tmp isa OpenBB.LinearConstraintSet
        @assert OpenBB.get_dependency(tmp) == [[1, 2]]
    end
end

for o in obj
    tmp = copy(deepcopy(o))
    tmp = sparse(tmp)
    @assert OpenBB.get_numVariables(tmp) == 4
    if tmp isa OpenBB.LinearObjective
        @assert OpenBB.get_dependency(tmp) == [1,2,3,4]
    elseif tmp isa OpenBB.QuadraticObjective
        @assert  OpenBB.get_dependency(tmp) == [1, 2, 3, 4]
    end
end


for v in var
    for c in cns
        for o in obj
            tmp = copy(deepcopy(OpenBB.Problem(varSet=v,cnsSet=c,objFun=o)))
            tmp = sparse(tmp)
            @assert OpenBB.get_numVariables(tmp) == 4
            @assert OpenBB.get_numDiscrete(tmp) == 2
            @assert OpenBB.get_numConstraints(tmp) == 1
            @assert OpenBB.get_variableBounds(tmp) == OpenBB.get_bounds(tmp.varSet)
            @assert OpenBB.get_discreteIndices(tmp) == OpenBB.get_discreteIndices(tmp.varSet)
            @assert OpenBB.get_sos1Groups(tmp) == OpenBB.get_sos1Groups(tmp.varSet)
            @assert OpenBB.get_pseudoCosts(tmp) == OpenBB.get_pseudoCosts(tmp.varSet)
            @assert OpenBB.get_constraintBounds(tmp) == OpenBB.get_bounds(tmp.cnsSet)
            @assert OpenBB.get_objective_dependency(tmp) == OpenBB.get_dependency(o)
            @assert OpenBB.get_constraints_dependency(tmp) == OpenBB.get_dependency(c)
        end
    end
end
