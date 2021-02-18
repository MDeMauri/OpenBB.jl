# @Author: Massimo De Mauri <massimo>
# @Date:   2021-02-12T13:11:59+01:00
# @Email:  massimo.demauri@protonmail.com
# @Filename: nonlinear_example.jl
# @Last modified by:   massimo
# @Last modified time: 2021-02-12T13:23:40+01:00
# @License: LGPL-3.0
# @Copyright: {{copyright}}

# load ipopt interface
OpenBB.load_subsolver_interface("IPOPT")


# create mixed integer non-linear problem (with linear/quadratic objective)
varSet = OpenBB.VariableSet(vals=zeros(4),loBs=-Infs(4),upBs=Infs(4),dscIndices=[1,2])

# # This get transformed into a generic convex function...
# objFun = OpenBB.QuadraticObjective(Q=Matrix{Float64}(I,4,4),L=zeros(4))

# create objective functions
(Q,L) = (Matrix{Float64}(I,4,4),zeros(4))
function evalObjVal(x)
    return .5*x'*Q*x + L'x
end
function evalObjGrd(x)
    return Q*x + L
end

function evalObjHes(x)
    return Q
end

grdSparsity = spzeros(Bool,4)
hesSparsity = sparse([1,2,3,4],[1,2,3,4],true,4,4)

objFun = OpenBB.ConvexObjective(evalVal=evalObjVal,evalGrd=evalObjGrd,evalHes=evalObjHes,
                                grdSparsity=grdSparsity,hesSparsity=hesSparsity,
                                typeGrd=Vector{Float64},typeHes=Matrix{Float64})


# create constraint functions
function evalCnsVal(x)
    return [ x[1] + x[2],
             x[3]^2 + x[4]^2,
             x[3] - 10.0*x[1],
             x[4] -  5.0*x[2] ]
end
cnsLoBs = [1.0, -Inf, 0.0, 0.0]
cnsUpBs = [1.0,  1e2, Inf, Inf]

function evalCnsJcb(x)
    return Matrix([    1.0         1.0     0.0         0.0         ;
                       0.0         0.0     2.0*x[3]    2.0*x[4]    ;
                       -10.0       0.0     1.0         0.0         ;
                       0.0         -5.0    0.0         1.0         ])
end
jcbSparsity = sparse(Matrix([   true      true    false   false   ;
                                false     false   true    true    ;
                                true      false   true    false   ;
                                false     true    false   true    ]))

evalCnsHes = Vector{Function}(undef,4)
evalCnsHes = x->[  zeros(4,4),
                zeros(4,4),
                vcat(zeros(2,4),hcat(zeros(2,2),Matrix(2.0I,2,2))),
                zeros(4,4)]

hesSparsity = [sparse(falses(4,4)),
              sparse(falses(4,4)),
              sparse(vcat(falses(2,4),hcat(falses(2,2),Matrix(true*I,2,2)))),
              sparse(falses(4,4))]

cnsSet = OpenBB.ConvexConstraintSet(evalVal=evalCnsVal,evalJcb=evalCnsJcb,evalHes=evalCnsHes,
                                       jcbSparsity=jcbSparsity,hesSparsity=hesSparsity,
                                       typeJcb=Matrix{Float64},typeHes=Matrix{Float64},
                                       loBs=cnsLoBs,upBs=cnsUpBs)

problem = OpenBB.Problem(varSet=varSet,cnsSet=cnsSet,objFun=objFun)
workspace = OpenBB.setup(problem,OpenBB.BBsettings(verbose=true,subsolverSettings=OpenBB.IPOPTsettings()))
OpenBB.solve!(workspace)
