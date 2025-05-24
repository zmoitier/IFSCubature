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
        @assert all(all(S.ρ < 1) for S in ifs) "The matrices `A` must be contractions."
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

function fix_points(sas::SelfAffineSet{D,T,N}) where {D,T,N}
    return fix_point.(sas.ifs)
end

function diameter(sas::SelfAffineSet{D,T,N}) where {D,T,N}
    return 2 * sas.bounding_ball.radius
end

function smallest_radius(
    z::AbstractVector, ifs::Vector{AffineMap{D,T,N}}, p::Real=2
) where {D,T,N}
    if p ≈ 2
        return maximum(norm(S(z) - z, p) / (1 - S.ρ) for S in ifs)
    end
    return maximum(norm(S(z) - z, p) / (1 - opnorm(S.A, p)) for S in ifs)
end

"""Return the smallest bounding ball."""
function bounding_ball(ifs::Vector{AffineMap{D,T,N}}; k::Int=1) where {D,T,N}
    ams = deepcopy(ifs)
    ins = 1 ./ (1 .- [opnorm(S.A, 2) for S in ams])
    z = MVector{D,T}(sum(fix_point.(ifs)) ./ length(ifs))

    result = Optim.optimize(c -> _smallest_radius(c, ams, ins, 2), z)
    z, r = Optim.minimizer(result), Optim.minimum(result)

    for _ in 2:k
        ams = _combine(ifs, ams)
        ins = 1 ./ (1 .- [opnorm(S.A, 2) for S in ams])

        result = Optim.optimize(c -> _smallest_radius(c, ams, ins, 2), z)
        z, r = Optim.minimizer(result), Optim.minimum(result)
    end

    return hyper_ball(z, r)
end

"""Return the smallest bounding box."""
function bounding_box(ifs::Vector{AffineMap{D,T,N}}; k::Int=1) where {D,T,N}
    op_norm = [opnorm(S.A, Inf) for S in ifs]
    for (S, n) in zip(ifs, op_norm)
        @assert !(isapprox(n, 1) || (n > 1)) "Affine map `$S` is not contracting for ∞-norm."
    end

    ams = deepcopy(ifs)
    ins = 1 ./ (1 .- op_norm)
    z = MVector{D,T}(sum(fix_point.(ifs)) ./ length(ifs))

    result = Optim.optimize(c -> _smallest_radius(c, ams, ins, Inf), z)
    z, r = Optim.minimizer(result), Optim.minimum(result)

    for _ in 2:k
        ams = _combine(ifs, ams)
        ins = 1 ./ (1 .- [opnorm(S.A, Inf) for S in ams])

        result = Optim.optimize(c -> _smallest_radius(c, ams, ins, Inf), z)
        z, r = Optim.minimizer(result), Optim.minimum(result)
    end

    return hyper_box(z, Diagonal(fill(r, D)))
end

function _combine(
    ifs_a::Vector{AffineMap{D,T,N}}, ifs_b::Vector{AffineMap{D,T,N}}
) where {D,T,N}
    ifs_c = Vector{AffineMap{D,T,N}}()
    for (R, S) in Iterators.product(ifs_a, ifs_b)
        push!(ifs_c, R ∘ S)
    end
    return ifs_c
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

    return (similarity_dimension(η), similarity_dimension(ρ))
end

function similarity_dimension(ifs::Vector{AffineMap{D,T,N}}) where {D,T,N}
    return similarity_dimension([S.ρ for S in ifs])
end

function similarity_dimension(ρ::Vector)
    bracket = (-log(length(ρ))) ./ log.(extrema(ρ))
    return find_zero(d -> sum(ρ .^ d) .- 1, bracket, Brent())
end
