module IFSCubature

import Base.Iterators: flatten
import Base: length, ∘
import IterTools: partition, product

using BoundingSphere
using FastGaussQuadrature
using IterativeSolvers
using LinearAlgebra
using StaticArrays

include("affine-map.jl")
include("shapes.jl")
include("self-affine-set.jl")
include("points.jl")
include("lagrange_polynomials.jl")
include("quadrature.jl")
include("example.jl")

end
