# @Author: Massimo De Mauri <massimo>
# @Date:   2019-05-12T21:19:38+02:00
# @Email:  massimo.demauri@gmail.com
# @Filename: multiprocessing_utils.jl
# @Last modified by:   massimo
# @Last modified time: 2019-06-06T12:52:39+02:00
# @License: LGPL-3.0
# @Copyright: {{copyright}}

# pause the current process without giving up the control of it
function pause(duration::Float)::Nothing
    start_time = time()
    while time() - start_time < duration
    end
    return
end

# define variables on remote workers
function remote_define(remoteSymbol::Symbol,variable::Any,processId::Int)::Nothing
    fetch(remote_do(OpenBB.eval,processId,Expr(:(=),remoteSymbol,:($variable))))
    return
end

function remote_define(remoteSymbols::Vector{Symbol},variables::Vector{T},processId::Int)::Nothing where T <: Any
    @assert length(remoteSymbols) == length(variables)
    for k in 1:length(remoteSymbols)
        remote_define(remoteSymbols[k],variables[k],processId)
    end
    return
end

function remote_define(remoteSymbol::Symbol,variable::Any,processIds::Vector{Int})::Nothing
    @sync for k in 1:length(processIds)
        @async remote_define(remoteSymbol,variable,processIds[k])
    end
    return
end

function remote_define(remoteSymbols::Vector{Symbol},variables::Vector{T},processIds::Vector{Int})::Nothing where T <: Any
    @sync for k in 1:length(processIds)
        @async remote_define(remoteSymbols,variables,processIds[k]::Int)
    end
    return
end
