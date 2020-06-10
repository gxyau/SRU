#===============================================================================
    Set of Range Uniqueness package. Main functionalities:
    - vandermonde generates a Vandermonde matrix using the base scalars (base)
      and dimension (dim)
    - Xmatrix generates X(S,π) for full rank linear maps/rank r linear maps
      where S is the proposed set for SRU and π is the permutation
===============================================================================#
module SRU

export vandermonde, Xmatrix

using Combinatorics, LinearAlgebra, Random

#===============================================================================
    column.(elm, dim) takes in a real number elm and
    an integer dim and generate (1, ..., elm^(dim - 1))
===============================================================================#
function column(elm::Real, dim::Integer)
    if dim <= 0
        error("Dimension has to be positive integer!")
    end

    return elm.^(1:dim)
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
        X[i,:] = [S[i] S̄[i]]
    end
    return X
end

end # SRU module
