# @Author: Massimo De Mauri <massimo>
# @Date:   2019-10-16T22:21:10+02:00
# @Email:  massimo.demauri@gmail.com
# @Filename: BBupdatesRegister.jl
# @Last modified by:   massimo
# @Last modified time: 2019-10-16T22:35:20+02:00
# @License: LGPL-3.0
# @Copyright: {{copyright}}

abstract type AbstractRegister end

struct BBupdatesRegister <: AbstractRegister
    changeList::Array{Function,1}
    pseudoCosts::Array{Float64,1}
end


# empty constructor
function BBupdatesRegister()::BBupdatesRegister
    return BBupdatesRegister(Function[],Float64[])
end
