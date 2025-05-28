module IFSCubature

import Base: length

using DataStructures: BinaryHeap
using FastGaussQuadrature: gausslegendre, gausslobatto
using IterTools: partition, product
using IterativeSolvers: powm!
using LinearAlgebra
using Roots: Brent, find_zero
using StaticArrays

const Float = float(Int)

include("affine-map.jl")
include("shapes.jl")
include("self-affine-set.jl")
include("points.jl")
include("lagrange-polynomials.jl")
include("cubature.jl")

include("example-1.jl")
include("example-2.jl")
include("example-3.jl")
include("example-n.jl")

end
