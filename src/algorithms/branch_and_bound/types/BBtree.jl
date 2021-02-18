# @Author: Massimo De Mauri <massimo>
# @Date:   2020-04-12T18:57:54+02:00
# @Email:  massimo.demauri@gmail.com
# @Filename: BBtree.jl
# @Last modified by:   massimo
# @Last modified time: 2021-02-14T17:22:33+01:00
# @License: LGPL-3.0
# @Copyright: {{copyright}}

const BBnodePool = Vector{Tuple{BBnode,Float}}

mutable struct BBtree
    active::BBnodePool
    solutions::BBnodePool
    suboptimals::BBnodePool
    infeasibles::BBnodePool
    blacklisted::BBnodePool
end

# empty constructor
function BBtree()::BBtree
    return BBtree(BBnodePool(),BBnodePool(),BBnodePool(),BBnodePool(),BBnodePool())
end

# empty the tree
function Base.empty!(tree::BBtree)::Nothing
    empty!(tree.active)
    empty!(tree.solutions)
    empty!(tree.suboptimals)
    empty!(tree.infeasibles)
    empty!(tree.blacklisted)
    return
end


function Base.sort!(pool::BBnodePool)::Nothing
    permute!(pool,sortperm(getindex.(pool,2),alg=MergeSort))
    return
end



# filtering function
function Base.filter!(selector::Function,tree::BBtree)::Nothing
    filter!(tuple->selector(tuple[1]),tree.active)
    filter!(tuple->selector(tuple[1]),tree.solutions)
    filter!(tuple->selector(tuple[1]),tree.infeasibles)
    filter!(tuple->selector(tuple[1]),tree.suboptimals)
    filter!(tuple->selector(tuple[1]),tree.blacklisted)
    return
end

# applies the given function to all the nodes in the tree
function applytoall!(fun::Function,tree::BBtree)::Nothing
    @. fun(getindex(tree.active,1))
    @. fun(getindex(tree.solutions,1))
    @. fun(getindex(tree.infeasibles,1))
    @. fun(getindex(tree.suboptimals,1))
    @. fun(getindex(tree.blacklisted,1))
    return
end
