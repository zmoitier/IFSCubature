"""Return the cantor dust for dimension ≥ 2."""
function cantor_dust(
    contraction_factor::T,
    fix_points_1d::AbstractVector{T},
    space_dim::Int,
    name::Union{String,Nothing}=nothing,
) where {T}
    @assert 0 ≤ contraction_factor ≤ 1 "the contraction factor must be between [0, 1)."
    @assert length(fix_points_1d) ≥ 1 "There must be at least 1 fix points."
    @assert space_dim ≥ 2 "the dimension must be greater or equal than 2."

    L = length(fix_points_1d)
    D = space_dim
    N = D * D

    fix_pts = sort(fix_points_1d)

    _cantor_open_set_condition(fill(contraction_factor, length(fix_points_1d)), fix_pts)

    ifs::Vector{AffineMap{D,T,N}} = []
    for c in product([fix_pts for _ in 1:space_dim]...)
        push!(ifs, contractive_similarity(contraction_factor, [v for v in c]))
    end

    ball = hyper_ball(
        fill((fix_pts[end] + fix_pts[1]) / 2, D), √D * (fix_pts[end] - fix_pts[1]) / 2
    )
    box = hyper_box(
        fill((fix_pts[end] + fix_pts[1]) / 2, D),
        Matrix(Diagonal(fill((fix_pts[end] - fix_pts[1]) / 2, D))),
    )

    if isnothing(name)
        name = "$(space_dim)d-cantor-dust"
    end

    return SelfAffineSet(ifs, fill(T(1//L^D), L^D), ball, box, name)
end
