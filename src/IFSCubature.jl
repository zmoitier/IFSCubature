module IFSCubature

using LinearAlgebra
using StaticArrays

import Base: length
import FastGaussQuadrature: gausslegendre, gausslobatto
import IterativeSolvers: powm!
import IterTools: partition, product
import Roots: Brent, find_zero

using Optim: Optim

include("affine-map.jl")
include("shapes.jl")
include("self-affine-set.jl")
include("points.jl")
include("lagrange_polynomials.jl")
include("quadrature.jl")

include("example-1.jl")
include("example-2.jl")
include("example-3.jl")
include("example-n.jl")

end
