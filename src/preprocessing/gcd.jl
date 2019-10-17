# @Author: Wim Van Roy
# @Date:   Mi Sep 25 14:26:39 CEST 2019
# @Email:  wim.vanroy@kuleuven.be
# @Filename: gcd.jl
# @Last modified by:   massimo
# @Last modified time: Mi Sep 25 14:26:54 CEST 2019
# @License: LGPL-3.0
# @Copyright: {{copyright}}


# Perform Greates Common Divider Update
function preprocess_gcd!( A::SparseMatrixCSC{Float64,Int},
                            cnsLoBs::Array{Float64,1},cnsUpBs::Array{Float64,1},
                            varLoBs::Array{Float64,1},varUpBs::Array{Float64,1},
                            dscIndices::Array{Int64,1})::Bool
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
            if !any(x->isinteger(x), multipliers)
                fractions = filter( x -> x != 0, unique(rem.(multipliers, 1)))
                commonMultiple = lcm.(Int.(floor.(inv.(p)+eps())))
                multipliers = commonMultiple * multipliers
                if any(x->isinteger(x), multipliers)
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
