# @Author: Massimo De Mauri <massimo>
# @Date:   2021-01-28T12:56:58+01:00
# @Email:  massimo.demauri@protonmail.com
# @Filename: sorting.jl
# @Last modified by:   massimo
# @Last modified time: 2021-01-28T14:28:21+01:00
# @License: LGPL-3.0
# @Copyright: {{copyright}}

function collect_score_and_sort!(V::Vector{T},scoreFun::Function;algorithm::Base.Sort.Algorithm=QuickSort,reverse::Bool=false)::Nothing where T
    permute!(V,sortperm(scoreFun.(V),rev=reverse))
    return
end
