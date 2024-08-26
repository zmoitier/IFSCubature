module IFSCubature

import Base.Iterators: flatten
import Base: length, ∘
import IterTools: partition, product
import Roots: Brent, find_zero

using FastGaussQuadrature
using IterativeSolvers
using LinearAlgebra
using StaticArrays

include("shapes.jl")

include("points.jl")

include("transformation.jl")

include("attractor.jl")

include("lagrange_polynomials.jl")

include("quadrature.jl")

include("example.jl")

end
