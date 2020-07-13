#===============================================================================
    Set of Range Uniqueness package. Main functionalities:
    - vandermonde generates a Vandermonde matrix using the base scalars (base)
      and dimension (dim)
    - Xmatrix generates X(S,π) for full rank linear maps/rank r linear maps
      where S is the proposed set for SRU and π is the permutation
===============================================================================#
module SRU

export vandermonde, Xmatrix, permutations_cycles_form

using Combinatorics, LinearAlgebra, Random

#===============================================================================
    column.(elm, dim) takes in a real number elm and
    an integer dim and generate (1, ..., elm^(dim - 1))
===============================================================================#
function column(elm::Real, dim::Integer)
    if dim <= 0
        error("Dimension has to be positive integer!")
    end

    return elm.^(0:dim-1)
end

#===============================================================================
    vandermonde(base, dim) takes in an array of numbers base, and an
    integer dim, and generate Vandermonde matrix with the specified
    base values.
===============================================================================#
function vandermonde(base, dim)
    if length(base) != dim
        throw("Number of base values provided is different than the dimension.")
    end

    return reduce(hcat, column.(base, dim))
end

#===============================================================================
    Xmatrix takes in a set of vectors, S, and a permutation, π, and generates
    the X(S,π) matrix for rank r linear maps.
===============================================================================#
function Xmatrix(S::Array{Array{Int64,1},1}, π::Array{Int64,1})
    k = length(π)
    n = length(S[1])
    S̄ = S[π]
    X = Array{Int64, 2}(undef, k, 2n)
    for i in 1:k
        X[i,:] = [S[i] -S̄[i]]
    end
    return X
end

function permutations_cycles_form(π::Array{Int64,1})
    cycleform = Array{Array{Int64,1}, 1}()
    n         = length(π)
    visited   = Bool[false for i=1:n]
    for i in 1:n
        if !visited[i]
            # Tortoise and hare cycle detection algorithm
            cycle         = Int64[]
            fast          = π[i]
            slow          = i
            push!(cycle, slow)
            while slow != fast
                slow          = π[slow]
                fast          = π[π[fast]]
                push!(cycle, slow)
            end
            push!(cycleform, cycle)
            # Update seen value
            for i in cycle
                visited[i] = true
            end
        end
    end
    return cycleform
end

end # SRU module
