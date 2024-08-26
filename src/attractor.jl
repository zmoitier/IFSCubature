struct ContractiveSimilarity{D,T,N}
    A::SMatrix{D,D,T,N}
    b::SVector{D,T}
    factor::T
end

function ContractiveSimilarity(
    factor::Real, matrix::Matrix{<:Real}, fix_point::Vector{<:Real}
)
    @assert 0 ≤ factor < 1 "0 ≤ factor=$factor < 1."
    @assert isapprox(matrix' * matrix, I) "$matrix must be an orthogonal matrix."

    D = length(fix_point)
    T = typeof(factor)
    N = D * D

    A = factor .* SMatrix{D,D,T,N}(matrix)
    b = (I - A) * SVector{D,T}(fix_point)

    return ContractiveSimilarity{D,T,N}(A, b, factor)
end

function ContractiveSimilarity(factor::Real, fix_point::Vector{<:Real})
    @assert 0 ≤ factor < 1 "0 ≤ factor=$factor < 1."

    D = length(fix_point)
    T = typeof(factor)
    N = D * D

    A = SMatrix{D,D,T,N}(factor .* I)
    b = (1 - factor) * SVector{D,T}(fix_point)

    return ContractiveSimilarity{D,T,N}(A, b, factor)
end

function (S::ContractiveSimilarity{1,T,1})(x::T)::T where {T<:Real}
    return S.A[1, 1] * x + S.b[1]
end

function (S::ContractiveSimilarity{D,T,N})(
    v::SVector{D,T}
)::SVector{D,T} where {D,T<:Real,N}
    return S.A * v + S.b
end

function fix_point(S::ContractiveSimilarity{D,T,N})::SVector{D,T} where {D,T,N}
    return (I - S.A) \ S.b
end

function (S::ContractiveSimilarity{D,T,N})(
    ball::HyperBall{D,T}
)::HyperBall{D,T} where {D,T,N}
    return HyperBall{D,T}(S.A * ball.center .+ S.b, S.factor * ball.radius)
end

function (S::ContractiveSimilarity{D,T,N})(
    box::HyperBox{D,T,N}
)::HyperBox{D,T,N} where {D,T,N}
    return HyperBox{D,T,N}(S.A * ball.center .+ S.b, S.A * box.paxis)
end

struct Attractor{D,T,N}
    name::String
    ifs::Vector{ContractiveSimilarity{D,T,N}} # Iterated Function System
    hausdorff_dim::Float64
    bounding_ball::HyperBall{D,T}
    bounding_box::HyperBox{D,T,N}
end
