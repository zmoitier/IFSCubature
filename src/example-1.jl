"""Return the Cantor set."""
function cantor_set(
    contraction_factors::AbstractVector,
    fix_points::AbstractVector,
    measure::AbstractVector,
    name::String="1d-cantor-set",
)
    L = length(fix_points)

    idx_sort = sortperm(fix_points)
    ρ = contraction_factors[idx_sort]
    c = fix_points[idx_sort]

    @assert L ≥ 1 "There must be at least 1 fix points."
    @assert L == length(ρ) "contraction_factors and fix_points must have the same length."
    @assert all(0 .≤ ρ .< 1) "the contraction factors must be between [0, 1)."

    _cantor_open_set_condition(ρ, c)

    ifs = [contractive_similarity(ρ, [c]) for (ρ, c) in zip(ρ, c)]

    b_ball = hyper_ball((c[end] + c[1]) / 2, (c[end] - c[1]) / 2)
    b_box = hyper_box((c[end] + c[1]) / 2, (c[end] - c[1]) / 2)

    return SelfAffineSet(ifs, measure, b_ball, b_box, name)
end

"""Return the Cantor set."""
function cantor_set(
    contraction_factor::Real, fix_points::AbstractVector, name::String="1d-cantor-set"
)
    T = typeof(contraction_factor)
    L = length(fix_points)
    return cantor_set(
        fill(contraction_factor, length(fix_points)), fix_points, fill(T(1//L), L), name
    )
end

"""Return the segment [a, b]."""
function segment(a::T, b::T) where {T}
    return cantor_set(T(1//2), [a, b], "1d-segment")
end

"""Warn if the open set condition is not satisfy."""
function _cantor_open_set_condition(ρ::AbstractVector{<:Real}, c::AbstractVector{<:Real})
    a, b = c[1], c[end]
    if !all(
        ρᵢ * (b - cᵢ) + cᵢ <= ρⱼ * (a - cⱼ) + cⱼ for
        ((ρᵢ, ρⱼ), (cᵢ, cⱼ)) in zip(partition(ρ, 2, 1), partition(c, 2, 1))
    )
        @warn "open set condition not satisfy."
    end

    return nothing
end
