# @Author: Massimo De Mauri <massimo>
# @Date:   2020-12-07Float18:20:42+01:00
# @Email:  massimo.demauri@gmail.com
# @Filename: BBblackList.jl
# @Last modified by:   massimo
# @Last modified time: 2020-12-17T15:43:34+01:00
# @License: LGPL-3.0
# @Copyright: {{copyright}}

mutable struct BBblackList
    assLen::Int
    assignments::Vector{Vector{Float}}
    splitFunc::Function
    function BBblackList(assLen::Int)::BBblackList
        return new(assLen,Vector{Vector{Real}}(),x->Int[])
    end
end


function Base.insert!(blackList::BBblackList,newAssignment::Vector{Float})::Nothing
    @assert length(newAssignment)==blackList.assLen
    if !lookup(blackList,newAssignment)
        push!(blackList.assignments,newAssignment)
        blackList.splitFunc = findsplit(blackList.assignments)
    end
    return
end

function Base.insert!(blackList::BBblackList,newAssignments::Vector{Vector{Float}})::Nothing
    @assert all(@. length(newAssignments)==blackList.assLen)
    newAssignments = unique(newAssignments)
    newAssignments = filter(x->!lookup(blackList,x),newAssignments)
    if !isempty(newAssignments)
        append!(blackList.assignments,newAssignments)
        blackList.splitFunc = findsplit(blackList.assignments)
    end
    return
end

function lookup(blackList::BBblackList,assignment::Vector{Float})::Bool
    index = blackList.splitFunc(assignment)
    if isempty(index)
        return false
    else
        toCompare = blackList.assignments[index[1]]
        if all(@. assignment == toCompare)
            return true
        else
            return false
        end
    end
end


function findsplit(assignments::Union{SubArray{Vector{Float},1},Vector{Vector{Float}}})::Function

    if length(assignments)==1
        return function locate0(assignment::Vector{Float})::Vector{Int}
                   return [1]
               end
    end

    numel=length(assignments[1])
    bestScore = -Inf
    bestIndex = 0
    bestSplit = (Int[],Int[])
    bestThreshold = Inf

    for k in 1:numel
        values = @. getindex(assignments,k)
        median_ = median(values)
        threshold = median_ + 0.1*sign(mean(values)-median_) # avoids problems with, for instace, [0,1,1]
        question = x->x>=threshold
        upGroup = findall(question,values)
        dnGroup = findall(!question,values)
        score = length(upGroup)*length(dnGroup)
        if score == 0.25*numel^2
            bestScore = 0.25*numel^2
            bestIndex = k
            bestThreshold = threshold
            bestSplit = (upGroup,dnGroup)
            break
        elseif score > bestScore
            bestScore = score
            bestIndex = k
            bestThreshold = threshold
            bestSplit = (upGroup,dnGroup)
        end
    end

    if any(isempty.(bestSplit))
        println(hcat(assignments...))
        println(bestIndex)
        println(bestThreshold)
        println(bestScore)
        error("BBblacklist: Splitting Threshold not Found")
    elseif length(bestSplit[1])==1 && length(bestSplit[2])<=1
        function locate1(assignment::Vector{Float})::Union{Vector{Int},SubArray{Int,1}}
            if assignment[bestIndex] >= bestThreshold
                return bestSplit[1]
            else
                return bestSplit[2]
            end
        end
        return locate1
    elseif length(bestSplit[1])==1
        nested = findsplit(view(assignments,bestSplit[2]))
        function locate2(assignment::Vector{Float})::Union{Vector{Int},SubArray{Int,1}}
            if assignment[bestIndex] >= bestThreshold
                return bestSplit[1]
            else
                return view(bestSplit[2],nested(assignment))
            end
        end
        return locate2
    elseif length(bestSplit[2])<=1
         nested = findsplit(view(assignments,bestSplit[1]))
            function locate3(assignment::Vector{Float})::Union{Vector{Int},SubArray{Int,1}}
                if assignment[bestIndex] >= bestThreshold
                    return view(bestSplit[1],nested(assignment))
                else
                    return bestSplit[2]
                end
            end
            return locate3
    else
        nested1 = findsplit(view(assignments,bestSplit[1]))
        nested2 = findsplit(view(assignments,bestSplit[2]))
        function locate4(assignment::Vector{Float})::Union{Vector{Int},SubArray{Int,1}}
            if assignment[bestIndex] >= bestThreshold
                return view(bestSplit[1],nested1(assignment))
            else
                return view(bestSplit[2],nested2(assignment))
            end
        end
        return locate4
    end
end




# test
# Len=100
# numData = 100
# data = [floor.(10*rand(Len)) for k in 1:numData]
# blackList = BBblackList(Len)
# insert!(blackList,data)
#
#
# for d in data
#     @info lookup(blackList,d)
# end
# println("-------------------------")
# for k in 1:2*numData
#     @time lookup(blackList,floor.(10*rand(Len)))
# end
