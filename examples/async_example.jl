# @Author: Massimo De Mauri <massimo>
# @Date:   2019-02-18T16:17:27+01:00
# @Email:  massimo.demauri@gmail.com
# @Filename: async_example.jl
# @Last modified by:   massimo
# @Last modified time: 2019-05-06T18:49:13+02:00
# @License: apache 2.0
# @Copyright: {{copyright}}

# this would show how to make OpenBB run asynchronously
using OpenBB
using SparseArrays
using LinearAlgebra



# I need to find an harder problem
problem = OpenBB.Problem(objFun=OpenBB.QuadraticObj(Q=Matrix(2.0I,8,8,),L=[-.5,2.,2.,2.,2.,2.,2.,2.]),
                         cnsSet=OpenBB.LinearCns(A = vcat(hcat(ones(1,4),zeros(1,4)),hcat(zeros(1,4),ones(1,4))),
                                                 loBs=[1.,1],upBs=[1.,1]),
                         varSet=OpenBB.VariableSet(loBs=[-5.;-Infs(3);-5;-Infs(3)],upBs=[-5.;Infs(3);5.;Infs(3)],val=zeros(8),dscIndices=[1,2,3,4,5,6,7,8]))

# setting flag to true stops the BBalgorithm
flag = Base.RefValue(false)
function stopBB(workspace::OpenBB.BBworkspace)
    return flag[]
end


# setup the BB process
workspace = OpenBB.setup(problem,OpenBB.BBsettings(verbose=true,custom_stopping_rule=stopBB),OpenBB.OSQPsettings())


# start the BB a first time and then stop it
println("the BB process was started...")
println("-the status of the BB process is: ",workspace.status.description)
task = @async OpenBB.solve!(workspace)
sleep(0.001)
flag[] = true
result1 = fetch(task)
workspace1 = deepcopy(workspace)
println("...the BB process was stopped.")
println("-the status of the BB process is: ",workspace.status.description,"\n")

# here we can change the constraints and bounds



# restart the BB process
println("The BB process was started again...")
flag[] = false
task = @async OpenBB.solve!(workspace)
println("-the status of the BB process is: ",workspace.status.description)
result2 = fetch(task)
workspace2 = deepcopy(workspace)
println("...the BB process is over")
println("-the status of the BB process is: ",workspace.status.description)
