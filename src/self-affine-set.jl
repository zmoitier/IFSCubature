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
        @assert all(all(opnorm(S.A) < 1) for S in ifs) "The matrices `A` must be contrations."
        @assert all(0 .≤ measure .< 1) "The measure weights must be in the interval (0, 1)."
        @assert sum(measure) ≈ 1 "The measure weights must sum to 1."

        return new{D,T,N}(ifs, measure, bounding_ball, bounding_box, name)
    end
end

function contractive_similarity(
    factor::Real, matrix::AbstractMatrix, fix_point::AbstractVector
)
    @assert 0 ≤ factor < 1 "0 ≤ factor=$factor < 1."
    @assert isapprox(matrix' * matrix, I) "$matrix must be an orthogonal matrix."

    A = factor .* matrix
    b = (I - A) * fix_point

    return affine_map(A, b)
end

function contractive_similarity(factor::Real, fix_point::AbstractVector)
    D = length(fix_point)
    return contractive_similarity(factor, Matrix(I, D, D), fix_point)
end

function fix_point(f::AffineMap{D,T,N}) where {D,T,N}
    return (I - f.A) \ f.b
end

function smallest_radius(
    z::AbstractVector, ifs::Vector{AffineMap{D,T,N}}, p::Real=2
) where {D,T,N}
    return maximum(norm(S(z) - z, p) / (1 - opnorm(S.A, p)) for S in ifs)
end

"""Return the smallest bounding ball."""
function bounding_ball(ifs::Vector{AffineMap{D,T,N}}) where {D,T,N}
    fs = 1 ./ (1 .- [opnorm(S.A) for S in ifs])

    x0 = MVector{D,T}(sum(fix_point.(ifs)) ./ length(ifs))
    res = Optim.optimize(z -> _smallest_radius(z, ifs, fs, 2), x0)
    return hyper_ball(Optim.minimizer(res), Optim.minimum(res))
end

"""Return the smallest bounding box."""
function bounding_box(ifs::Vector{AffineMap{D,T,N}}) where {D,T,N}
    opnorm_inf = [opnorm(S.A, Inf) for S in ifs]
    for (S, n) in zip(ifs, opnorm_inf)
        @assert !(isapprox(n, 1) || (n > 1)) "Affine map `$S` is not contracting for ∞-norm."
    end

    fs = 1 ./ (1 .- opnorm_inf)
    x0 = MVector{D,T}(sum(fix_point.(ifs)) ./ length(ifs))
    res = Optim.optimize(z -> _smallest_radius(z, ifs, fs, Inf), x0)

    return hyper_box(Optim.minimizer(res), Diagonal(fill(Optim.minimum(res), D)))
end

function _smallest_radius(
    z::Union{SVector{D,T},MVector{D,T}},
    ifs::Vector{AffineMap{D,T,N}},
    fs::Vector{T},
    p::Real,
) where {D,T,N}
    return maximum(f * norm(S(z) - z, p) for (S, f) in zip(ifs, fs))
end

function dimension(ifs::Vector{AffineMap{D,T,N}}) where {D,T,N}
    η, ρ = zeros(T, length(ifs)), zeros(T, length(ifs))
    for (i, S) in enumerate(ifs)
        σ = svdvals(S.A)
        η[i], ρ[i] = σ[end], σ[1]
    end

    return (_similarity_dimension(η), _similarity_dimension(ρ))
end

function similarity_dimension(ifs::Vector{AffineMap{D,T,N}}) where {D,T,N}
    return _similarity_dimension([opnorm(S.A) for S in ifs])
end

function _similarity_dimension(ρ::Vector)
    bracket = (-log(length(ρ))) ./ log.(extrema(ρ))
    return find_zero(d -> sum(ρ .^ d) .- 1, bracket, Brent())
end