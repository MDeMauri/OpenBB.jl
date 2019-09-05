# @Author: Massimo De Mauri <massimo>
# @Date:   2019-03-11T17:03:23+01:00
# @Email:  massimo.demauri@gmail.com
# @Filename: linear_bounds_propagation.jl
# @Last modified by:   massimo
# @Last modified time: 2019-09-02T13:34:08+02:00
# @License: LGPL-3.0
# @Copyright: {{copyright}}


# single constraint linear bounds propagation
function bounds_propagation!(row::Int,
                            A::SparseMatrixCSC{Float64,Int},
                            cnsLoBs::Array{Float64,1},cnsUpBs::Array{Float64,1},
                            varLoBs::Array{Float64,1},varUpBs::Array{Float64,1},
                            dscIndices::Array{Int64,1})::Set{Int}
      cns = A[row,1:end]
      cnsLoB = cnsLoBs[row]
      cnsUpB = cnsUpBs[row]

      # use a set for the output
      updatedVars = Set{Int}()

      if cnsLoB == Inf && cnsUpB == Inf
          return updatedVars
      end

      # get the sparsity of cns
      indices,coeffs = findnz(cns)

      # distiguish between positive and negative coefficients
      posCoeffs = Array{Float64,1}(undef,length(indices))
      @. posCoeffs = coeffs*(coeffs>0)
      negCoeffs = Array{Float64,1}(undef,length(indices))
      @. negCoeffs = coeffs*(coeffs<0)

      # compute the max and the min value of the constraint components
      maxArray = Array{Float64,1}(undef,length(indices))
      @. maxArray = posCoeffs*varUpBs[indices] + negCoeffs*varLoBs[indices]
      minArray = Array{Float64,1}(undef,length(indices))
      @. minArray = posCoeffs*varLoBs[indices] + negCoeffs*varUpBs[indices]

      totalMaxArray = sum(maxArray)
      totalMinArray = sum(minArray)

      if totalMaxArray <= cnsUpB && totalMinArray >= cnsLoB
          # redundant constraint -> Remove?
          cnsLoBs[row] = -Inf
          cnsUpBs[row] = Inf
          return updatedVars
      elseif totalMinArray > cnsUpBs[row]
          throw(InfeasibleError("Infeasible bound detected"))
      elseif totalMaxArray < cnsLoBs[row]
          throw(InfeasibleError("Infeasible bound detected"))
      end


      # column index
      iteration = 0
      lastChanged = 1
      while true

            # select the variable to work on
            i = mod(iteration,length(indices))+1

            # compute new lower and upper bounds
            if coeffs[i] > 0
                  newLoB = (-(totalMaxArray - maxArray[i]) + cnsLoB)/coeffs[i]
                  newUpB = (-(totalMinArray - minArray[i]) + cnsUpB)/coeffs[i]
            else
                  newUpB = (-(totalMaxArray - maxArray[i]) + cnsLoB)/coeffs[i]
                  newLoB = (-(totalMinArray - minArray[i]) + cnsUpB)/coeffs[i]
            end

            if newLoB > newUpB
                  # Found infeasibility!
                  @info i, indices[i], varLoBs, varUpBs, newLoB, newUpB
                  throw(InfeasibleError("Infeasible bound detected"))
            end

            # perform changes
            if varLoBs[indices[i]] < newLoB
                  if indices[i] in dscIndices
                      newLoB = ceil(newLoB)
                  end

                  # change max and min arrays
                  if coeffs[i] > 0
                      totalMinArray = totalMinArray - minArray[i]
                      minArray[i] = coeffs[i]*newLoB
                      totalMinArray = totalMinArray + minArray[i]
                  else
                      totalMaxArray = totalMaxArray - maxArray[i]
                      maxArray[i] = coeffs[i]*newLoB
                      totalMaxArray = totalMaxArray + maxArray[i]
                  end
                  # update the bound
                  varLoBs[indices[i]] = newLoB
                  # remember that there was an update
                  lastChanged = i
                  # collect the updated variables
                  push!(updatedVars,indices[i])
            end

            if varUpBs[indices[i]] > newUpB
                  if indices[i] in dscIndices
                      newLoB = floor(newUpB)
                  end

                  # change max and min arrays
                  if coeffs[i] > 0
                      totalMaxArray = totalMaxArray - maxArray[i]
                      maxArray[i] = coeffs[i]*newUpB
                      totalMaxArray = totalMaxArray + maxArray[i]
                  else
                      totalMinArray = totalMinArray - minArray[i]
                      minArray[i] = coeffs[i]*newUpB
                      totalMinArray = totalMinArray + minArray[i]
                  end
                  # update the bound
                  varUpBs[indices[i]] = newUpB
                  # remember that there was an update
                  lastChanged = i
                  # collect the updated variables
                  push!(updatedVars,indices[i])
            end

            # If a variable bound changed, change cns bounds
            if lastChanged == i
                if totalMaxArray < cnsUpB
                    cnsUpB = totalMaxArray
                    cnsUpBs[row] = cnsUpB
                end
                if totalMinArray > cnsLoB
                    cnsLoB = totalMinArray
                    cnsLoBs[row] = cnsLoB
                end
            end

            # if an update took place:
            # mark the current variable as updated and restart the iteration
            if lastChanged - 1 == mod(i,length(indices))
                  return updatedVars
            else
                  iteration = iteration + 1
            end
      end
end


# multi-constraint linear bounds propagation
function bounds_propagation!(rowsToCheck::Set{Int},
                            A::SparseMatrixCSC{Float64,Int},
                            cnsLoBs::Array{Float64,1},cnsUpBs::Array{Float64,1},
                            varLoBs::Array{Float64,1},varUpBs::Array{Float64,1},
                            dscIndices::Array{Int64,1})::Set{Int}

      # use a set for the output
      updatedVars = Set{Int}()

      while length(rowsToCheck) > 0

            # pick the first variable to update
            row = pop!(rowsToCheck)

            # perform bound propagation and collect the updated variables
            newUpdatedVars = bounds_propagation!(row, A, cnsLoBs, cnsUpBs, varLoBs, varUpBs, dscIndices)

            if length(newUpdatedVars) > 0
                  # collect the new rows to check
                  newRowsToCheck = unique(findnz(A[1:end,collect(newUpdatedVars)])[1])
                  deleteat!(newRowsToCheck,findfirst(x->x==row,newRowsToCheck)) # do not reinsert the current row
                  union!(rowsToCheck,newRowsToCheck)
                  # remember which variables were updated
                  union!(updatedVars,newUpdatedVars)
            end
      end

      return updatedVars
end

function bounds_propagation!(node::BBnode, A::SparseMatrixCSC{Float64,Int}, dscIndices::Array{Int64,1}, updatedVars::Array{Int64,1})::Set{Int}
      if 0 in updatedVars
            rows = Set(1:size(A)[1])
      else
            rows = Set(unique(findnz(A[1:end,collect(updatedVars)])[1]))
      end

      return bounds_propagation!(rows,
                                 A, node.cnsLoBs, node.cnsUpBs,
                                 node.varLoBs, node.varUpBs, dscIndices)
end
