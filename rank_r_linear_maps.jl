#============================================================================
    On invertible matrices
============================================================================#
# Set current path
cd("/home/gxyau/Documents/github/SRU")

# Use package for linear algebra, custom, etc
using Revise
using Combinatorics, LinearAlgebra, Random, RowEchelon
include("./SRU.jl")
using .SRU

#============================================================================
    Computation for full ranked linear maps. Want to show that
    rank(X(S,π)) ≥ n+1 for all π ≠ id
============================================================================#
# Defining basic parameters
n   = 4 # Dimension of input vector space
r   = 3
k   = 2n-r+1
In  = Matrix(I,n,n) * 1 # Identity matrix
# Generates x1, ..., xk
for i in 1:n
    @eval $(Symbol(:x, i)) = In[:,$i]
end

for i in 1:(k-n)
    @eval $(Symbol(:x, n+i)) = (collect(1:n).^($i))
end

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
    X̂      = Xmatrix(S,perms[i])
    r̂      = rank(X̂)
    cycles = permutations_cycles_form(perms[i])
    if (r̂ <= 2n-r)
        println("The matrix is:")
        display(X̂)
        println("Permutation $(perms[i]) has rank $r̂")
        println("There are $(length(cycles)) cycles in this permutation")
        println(cycles)
        println("\n")
    end
end

#============================================================================
    Computation for full ranked linear maps. Want to show that
    rank(X(S,π)) ≥ n+1 for all π ≠ id
============================================================================#
n  = 5 # Dimension of input vector space
r  = 3 # Rank of the matrix
k  = 2n-r+1 # Number of vectors in S
V1 = Matrix(I,n,n) * 1
V2 = vandermonde(Int64[2,3,5,7,11],n)

# Generating x1, ... xk
for i in 1:n
    @eval $(Symbol(:x, i))   = V1[:,$i]
    @eval $(Symbol(:x, n+i)) = V2[:,$i]
end

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

# Number of permutations that violates the rank assumption
count        = 0
count_cycles = 0

# Compute ranks for X(S,π) for π ∈ Sym(k)
for i in 2:length(perms) # Can start from 1 but that's identity
    global count, count_cycles
    X̂ = Xmatrix(S,perms[i])
    r̂ = rank(X̂)
    if (r̂ <= 2n-2)
        count += 1
        println("The matrix is:")
        display(X̂)
        cycles = permutations_cycles_form(perms[i])
        println("Permutation $(perms[i]) has rank $r̂")
        println("There are $(length(cycles)) cycles in this permutation")
        println(cycles)
        println("\n")
        if (length(cycles) < n)
            count_cycles += 1
        end
    end
end

println("There are $count many permutations of not full rank")
println("There are $count_cycles many permutations of < $n cycles")
