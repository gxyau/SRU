#===================================================================
    Script to figure out whether or not a matrix is invertible
===================================================================#
# Set current path
cd("/home/gxyau/Documents/github/SRU")

# Use package for linear algebra, custom, etc
include("./SRU.jl")
using LinearAlgebra, Random, RowEchelon
using .SRU

# Setting seed
Random.seed!(1)

# Aliasing packages
const van = Vandermonde

#=
# Testing main functionalities of Vandermonde
points1 = van.column.(1:4, 4)
points2 = van.column.(5:8, 4)
matrix1 = van.columnstomatrix(points2)
matrix2 = van.columnstomatrix(points2)
=#

# Linear map from with matrix representation having m * n dimension
(mlinear,nlinear) = (1,2) # Dimensions of linear map/matrix

# Getting blocks

# Simulate permutations and generate determinants
nsim = 100
determinants = zeros(nsim)
for iteration in 1:nsim
    println("Iteration $iteration")
    A, B, C, D = generateblocks(nlinear)
    determinants[iteration] += reduceddeterminant(A, B, C, D)
    if determinants[iteration] == 0
        println([A B; C D])
    end
end
any(determinants .== 0)

A, B, C, D = generateblocks(nlinear)
[A B; C D]
reduceddeterminant(A, B, C, D)


#############################################################################
c = Int64[11, -16]
d = Int64[0,0]
c0 = 35
d0 = 1
s1, s2, s3, s4 = Int64[3, 2], Int64[-5, 1], Int64[11, 0], Int64[3, 14]
X1 = [hcat(Int64[1], transpose(s1)); hcat(Int64[1], transpose(s3))]
X2 = [hcat(Int64[1], transpose(s2)); hcat(Int64[1], transpose(s4))]
X3 = [hcat(Int64[1], transpose(s1)); hcat(Int64[1], transpose(s4))]
X4 = [hcat(Int64[1], transpose(s2)); hcat(Int64[1], transpose(s4))]
function f1(x)
    return d0 + c0 + (d+c) ⋅ x
end

function f2(x)
    return d0 - c0 + (d-c) ⋅ x
end

############################################################################
# Test
############################################################################
S = [s1,s2,s3,s4]
println(f1.(S))
println(f2.(S))

############################################################################
ctilde0 = 37
dtilde0 = 0
ctilde  = Int64[25; -8]
dtilde  = Int64[-1; 8]
f3(x) = dtilde0 + ctilde0 + (dtilde + ctilde) ⋅ x
f4(x) = dtilde0 - ctilde0 + (dtilde - ctilde) ⋅ x
Stilde = [s1, s2, s4]

println(f3.(Stilde))
println(f4.(Stilde))

############################################################################
s1, s2, s3, s4 = Float64[√2, 0, √2], Float64[0, √2, 0], Float64[-1, -1, 1], Float64[0, 0, 1]
X1 = [transpose(s1 + s2); transpose(s3 + s4)]
X2 = [transpose(s1 - s2); transpose(s3 - s4)]
c, d = Int64[1, -1, 0], Int64[1, -1, -2]
f(x) = (d+c) ⋅ x
g(x) = (d-c) ⋅ x
broadcast(f, [s1,s2,s3,s4])
broadcast(g, [s1,s2,s3,s4])
