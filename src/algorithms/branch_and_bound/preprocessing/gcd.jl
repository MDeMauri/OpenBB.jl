# @Author: Wim Van Roy
# @Date:   Mi Sep 25 14:26:39 CEST 2019
# @Email:  wim.vanroy@kuleuven.be
# @Filename: gcd.jl
# @Last modified by:   massimo
# @Last modified time: 2020-11-06T14:46:24+01:00
# @License: LGPL-3.0
# @Copyright: {{copyright}}


# Perform Greatest Common Divider Update
function preprocess_gcd!(A::SpMatrix{Float},
                         cnsLoBs::Vector{Float},cnsUpBs::Vector{Float},
                         varLoBs::Vector{Float},varUpBs::Vector{Float},
                         dscIndices::Vector{Int64})::Bool
    apply::Bool=false

    # Search for suitable rows
    rowsToCheck = unique(findnz(A[1:end,collect(dscIndices)])[1])

    # Go over suitable rows and find gcd
    while length(rowsToCheck) > 0
        # pick the first variable to update
        row = pop!(rowsToCheck)

        # Check if subset of dscIndices, if so calculate ggd & perform update
        items = findnz(A[row, 1:end])
        indices = unique(items[1])
        multipliers = unique(items[2])
        if issubset(indices, dscIndices)
            commonDivider = 1
            commonMultiple = 1
            apply = false

            # Make multipliers integer
            if any(x->(!isinteger(x)), multipliers)
                fractions = filter( x -> x != 0, unique(rem.(multipliers, 1)))
                inverse = inv.(fractions).+eps()
                commonMultiple = lcm(Int.(floor.(inverse)))
                multipliers = commonMultiple * multipliers
                if !any(x->(!isinteger(x)), multipliers)
                    apply = true
                end
            else
                apply = true
            end

            if apply
                commonDivider = gcd(Int.(multipliers))

                if commonDivider > 1
                    proportion = commonMultiple / commonDivider

                    # Update bounds
                    cnsLoBs[row] = ceil(cnsLoBs[row] * proportion - eps())
                    cnsUpBs[row] = floor(cnsUpBs[row] * proportion + eps())
                    A[row, 1:end] = A[row, 1:end] * proportion

                    if cnsLoBs[row] > cnsUpBs[row]
                          # Found infeasibility!
                          return false
                    end
                end
            end
        end
    end

    return true
end

function preprocess_gcd!(A::Matrix{Float},
                         cnsLoBs::Vector{Float},cnsUpBs::Vector{Float},
                         varLoBs::Vector{Float},varUpBs::Vector{Float},
                         dscIndices::Vector{Int64})::Bool
    return preprocess_gcd!(sparse(A),cnsLoBs,cnsUpBs,varLoBs,varUpBs,dscIndices)
end
