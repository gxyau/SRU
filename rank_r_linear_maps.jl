#===============================================================================
    On invertible matrices
===============================================================================#
# Set current path
cd("/home/gxyau/Documents/julia/research/")

# Use package for linear algebra, custom, etc
using Combinatorics, LinearAlgebra, Random
include("./SRU.jl")
using .SRU

#===============================================================================
    Computation for full ranked linear maps. Want to show that
    rank(X(S,π)) ≥ n+1 for all π ≠ id
===============================================================================#
# Defining basic parameters
n   = 4 # Dimension of input vector space
k   = n+1
In  = Matrix(I,n,n) * 1 # Identity matrix
# Generates x1, ..., xk
for i in 1:n
    @eval $(Symbol(:x, i)) = In[:,$i]
end
@eval $(Symbol(:x, k)) = collect(1:n)

# Putting the x's into S
S = Array{Array{Int64, 1}, 1}()
for i in 1:k
    @eval push!(S, $(Symbol(:x, i)))
end

#=
    This line generates all possible permutation on k symbols, one row
    for each permutation.
=#
perms = permutations(1:k, k) |> collect

# Compute ranks for X(S,π) for π ∈ Sym(k)
for i in 2:length(perms) # Can start from 1 but that's identity
    println("Permutation $(perms[i])")
    X̂ = Xmatrix(S,perms[i])
    r = rank(X̂)
    println(r)
end
