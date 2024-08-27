struct SelfAffineSet{D,T,N}
    ifs::Vector{AffineMap{D,T,N}} # Iterated Function System
    measure::Vector{T}
    bounding_ball::HyperBall{D,T}
    bounding_box::HyperBox{D,T,N}
    name::String

    function SelfAffineSet(
        ifs::Vector{AffineMap{D,T,N}},
        measure::Vector{T},
        bounding_ball::HyperBall{D,T},
        bounding_box::HyperBox{D,T,N},
        name::String,
    ) where {D,T,N}
        @assert length(ifs) == length(measure) "length(ifs) == length(measure)"
        @assert all(all(f.ρ .< 1) for f in ifs) "The matrices `A` must be contrations."
        @assert all(0 .< measure .< 1) "The measure weights must be in the interval (0, 1)."
        @assert sum(measure) ≈ 1 "The measure weights must sum to 1."

        return new{D,T,N}(ifs, measure, bounding_ball, bounding_box, name)
    end
end

function contractive_similarity(
    factor::Real, matrix::Matrix{<:Real}, fix_point::AbstractVector{<:Real}
)
    @assert 0 ≤ factor < 1 "0 ≤ factor=$factor < 1."
    @assert isapprox(matrix' * matrix, I) "$matrix must be an orthogonal matrix."

    A = factor .* matrix
    b = (I - A) * fix_point

    return affine_map(A, b)
end

function contractive_similarity(factor::Real, fix_point::AbstractVector{<:Real})
    @assert 0 ≤ factor < 1 "0 ≤ factor=$factor < 1."

    D = length(fix_point)
    A = factor .* Matrix(I, D, D)
    b = (I - A) * fix_point

    return affine_map(A, b)
end

function fix_point(f::AffineMap{D,T,N}) where {D,T,N}
    return (I - f.A) \ f.b
end

"""Return the bounding ball."""
function bounding_ball(ifs::Vector{AffineMap{D,T,N}}) where {D,T,N}
    ball_mean = _ball_mean(ifs)
    ball_fix_pts = _ball_fix_pts(ifs)

    if ball_mean.radius < ball_fix_pts.radius
        return ball_mean
    end

    return ball_fix_pts
end

function _ball_mean(ifs::Vector{AffineMap{D,T,N}}) where {D,T,N}
    s = sum(S.ρ for S in ifs)
    A = sum(S.ρ * S.A for S in ifs) ./ s
    b = sum(S.ρ * S.b for S in ifs) ./ s

    z = (I - A) \ b
    r = maximum([norm(S.A * z + S.b - z) / (1 - S.ρ) for S in ifs])

    return HyperBall{D,T}(z, r)
end

function _ball_fix_pts(ifs::Vector{AffineMap{D,T,N}}) where {D,T,N}
    c = [fix_point(S) for S in ifs]

    z, _ = boundingsphere(c)
    r = maximum([norm(S.A * z + S.b - z) / (1 - S.ρ) for S in ifs])

    return HyperBall{D,T}(z, r)
end
