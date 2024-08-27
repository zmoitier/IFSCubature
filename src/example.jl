#!#############################
#! Helper functions
#!#############################

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

#!#############################
#! 1-dimensional IFS
#!#############################

"""Return the segment [a, b]."""
function segment(a::T, b::T) where {T}
    @assert a ≤ b "must have a ≤ b."

    ρ = T(1//2)
    ifs = [contractive_similarity(ρ, [a]), contractive_similarity(ρ, [b])]

    b_ball = hyper_ball((a + b) / 2, (b - a) / 2)
    b_box = hyper_box((a + b) / 2, (b - a) / 2)

    return SelfAffineSet(ifs, fill(ρ, 2), b_ball, b_box, "1d-segment")
end

"""Return the Cantor set."""
function cantor_set(
    contraction_factors::AbstractVector{T}, fix_points::AbstractVector{T}
) where {T}
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

    return SelfAffineSet(ifs, fill(T(1//L), L), b_ball, b_box, "1d-cantor-set")
end

"""Return the Cantor set."""
function cantor_set(contraction_factor::T, fix_points::AbstractVector{T}) where {T}
    return cantor_set(fill(contraction_factor, length(fix_points)), fix_points)
end

#!#############################
#! 2-dimensional IFS
#!#############################

"""Return the square [a, b]x[a, b]."""
function square(a::T, b::T) where {T}
    return cantor_dust(T(1//2), [a, b], 2, "2d-square")
end

"""Return the triangle."""
function triangle()
    ifs = [contractive_similarity(0.5, Matrix(-I, 2, 2), fill(0.0, 2))]
    for c in [[1.0, 0.0], [-0.5, √3 / 2], [-0.5, -√3 / 2]]
        push!(ifs, contractive_similarity(0.5, c))
    end

    b_ball = hyper_ball(fill(0.0, 2), 1.0)
    b_box = hyper_box([0.25, 0], Matrix(Diagonal([0.75, √3 / 2])))

    return SelfAffineSet(ifs, fill(1 / 4, 4), b_ball, b_box, "2d-triangle")
end

"""Return the Sierpinski triangle."""
function sierpinski_triangle(contraction_factor::Union{Float64,Nothing}=nothing)
    if isnothing(contraction_factor)
        ρ = 0.5
    else
        ρ = contraction_factor
    end

    @assert 0 ≤ ρ < 1 "the contraction factor must be between [0, 1)."
    if ρ > 0.5
        @warn "open set condition not satisfy."
    end

    ifs = [
        contractive_similarity(ρ, c) for c in [[1.0, 0.0], [-0.5, √3 / 2], [-0.5, -√3 / 2]]
    ]

    b_ball = hyper_ball(fill(0.0, 2), 1.0)
    b_box = hyper_box([0.25, 0], Matrix(Diagonal([0.75, √3 / 2])))

    return SelfAffineSet(ifs, fill(1 / 3, 3), b_ball, b_box, "2d-sierpinski-triangle")
end

"""Return the 2d Vicsek fractal."""
function vicsek_2d(contraction_factor::Float64, angle::Union{Float64,Nothing}=nothing)
    @assert 0 ≤ contraction_factor ≤ 1 "the contraction factor must be between [0, 1)."

    name = "2d-vicsek"

    if isnothing(angle)
        Ri = Matrix(I, 2, 2)
    else
        Ri = matrix_rotation_2d(angle)
        name *= "-rot"
    end

    ifs::Vector{AffineMap{2,Float64,4}} = [
        contractive_similarity(contraction_factor, Ri, fill(0.0, 2))
    ]
    for c in product([-1, 1], [-1, 1])
        push!(ifs, contractive_similarity(contraction_factor, [v for v in c]))
    end

    b_ball = hyper_ball(fill(0.0, 2), √2)
    b_box = hyper_box(fill(0.0, 2), Matrix(Diagonal(fill(1.0, 2))))

    return SelfAffineSet(ifs, fill(1 / 5, 5), b_ball, b_box, name)
end

"""Return the Sierpinski carpet."""
function sierpinski_carpet()
    ρ = 1 / 3
    ifs = [
        contractive_similarity(ρ, c) for
        c in [[1, 0], [1, 1], [0, 1], [-1, 1], [-1, 0], [-1, -1], [0, -1], [1, -1]]
    ]

    b_ball = hyper_ball(fill(0.0, 2), √2)
    b_box = hyper_box(fill(0.0, 2), Matrix(Diagonal(fill(1.0, 2))))

    return SelfAffineSet(ifs, fill(1 / 8, 8), b_ball, b_box, "2d-sierpinski-carpet")
end

"""Return the Koch snowflake."""
function koch_snowflake()
    R = matrix_rotation_2d(1 / 6; implicit_pi=true)

    ifs = [contractive_similarity(1 / √3, R, fill(0.0, 2))]

    ρ = 1 / 3
    for (s, c) in sincospi.((0:5) .// 3)
        push!(ifs, contractive_similarity(ρ, [c, s]))
    end

    b_ball = hyper_ball(fill(0.0, 2), 1.0)
    b_box = hyper_box(fill(0.0, 2), Matrix(Diagonal([1.0, √3 / 2])))

    μ = [1 / 3, fill(1 / 9, 6)...]

    return SelfAffineSet(ifs, μ, b_ball, b_box, "2d-koch-snowflake")
end

"""Return the Gosper flowsnake."""
function gosper_flowsnake()
    ρ = 1 / √7
    R = matrix_rotation_2d(atan(√3 / 5))

    ifs = [contractive_similarity(ρ, R, fill(0.0, 2))]
    for (s, c) in sincospi.((0:5) .// 3)
        push!(ifs, contractive_similarity(ρ, R, [c, s]))
    end

    b_ball = hyper_ball(fill(0.0, 2), √3 / (√7 - 1))
    b_box = hyper_box(fill(0.0, 2), Matrix(Diagonal(fill(b_ball.radius, 2))))

    return SelfAffineSet(ifs, fill(1 / 7, 7), b_ball, b_box, "2d-gosper-flowsnake")
end

#!#############################
#! 3-dimensional IFS
#!#############################

"""Return the cube [a, b]x[a, b]x[a, b]."""
function cube(a::T, b::T) where {T}
    return cantor_dust(T(1//2), [a, b], 3, "3d-cube")
end

# """Return the tetrahedron."""
# function tetrahedron()
# end

"""Return the Sierpinski tetrahedron."""
function sierpinski_tetrahedron(contraction_factor::Union{Float64,Nothing}=nothing)
    if isnothing(contraction_factor)
        ρ = 0.5
    else
        ρ = contraction_factor
    end

    @assert 0 ≤ ρ ≤ 1 "the contraction factor must be between [0, 1)."

    if ρ > 0.5
        @warn "open set condition not satisfy."
    end

    ifs = [
        contractive_similarity(ρ, c) for
        c in [[1, 0, -√2 / 2], [-1, 0, -√2 / 2], [0, 1, √2 / 2], [0, -1, √2 / 2]]
    ]

    b_ball = hyper_ball(fill(0.0, 3), √(3 / 2))
    b_box = hyper_box(fill(0.0, 3), Matrix(Diagonal([1.0, 1.0, √2 / 2])))

    return SelfAffineSet(ifs, fill(1 / 4, 4), b_ball, b_box, "3d-sierpinski-tetrahedron")
end

"""Return Menger sponge."""
function menger_sponge()
    ρ = 1 / 3
    fix_points = [
        [1, 0, 1],
        [1, 1, 1],
        [0, 1, 1],
        [-1, 1, 1],
        [-1, 0, 1],
        [-1, -1, 1],
        [0, -1, 1],
        [1, -1, 1],
        #
        [1, 1, 0],
        [-1, 1, 0],
        [-1, -1, 0],
        [1, -1, 0],
        #
        [1, 0, -1],
        [1, 1, -1],
        [0, 1, -1],
        [-1, 1, -1],
        [-1, 0, -1],
        [-1, -1, -1],
        [0, -1, -1],
        [1, -1, -1],
    ]

    ifs = [contractive_similarity(ρ, c) for c in fix_points]

    b_ball = hyper_ball(fill(0.0, 3), √3)
    b_box = hyper_box(fill(0.0, 3), Matrix(Diagonal(fill(1.0, 3))))

    return SelfAffineSet(ifs, fill(1 / 20, 20), b_ball, b_box, "3d-menger-sponge")
end

"""Return the 3d Vicsek fractal."""
function vicsek_3d(contraction_factor::Float64, central_rotation::Bool=false)
    @assert 0 ≤ contraction_factor ≤ 1 "the contraction factor must be between [0, 1)."

    name = "3d-vicsek"
    if central_rotation
        Ry = matrix_rotation_3d(atan(√2), [0.0, 1.0, 0.0])
        Rz = matrix_rotation_3d(-0.25, [0.0, 0.0, 1.0]; implicit_pi=true)
        R = Rz * Ry
        name *= "-rot"
    else
        R = Matrix(I, 3, 3)
    end
    S0 = contractive_similarity(contraction_factor, R, fill(0.0, 3))

    ifs::Vector{AffineMap{3,Float64,9}} = [S0]
    for c in product([-1, 1], [-1, 1], [-1, 1])
        push!(ifs, contractive_similarity(contraction_factor, [v for v in c]))
    end

    b_ball = hyper_ball(fill(0.0, 3), √3)
    b_box = hyper_box(fill(0.0, 3), Matrix(Diagonal(fill(1.0, 3))))

    return SelfAffineSet(ifs, fill(1 / 9, 9), b_ball, b_box, name)
end

#!#############################
#! n-dimensional IFS
#!#############################

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

    b_ball = hyper_ball(
        fill((fix_pts[end] + fix_pts[1]) / 2, D), √D * (fix_pts[end] - fix_pts[1]) / 2
    )
    b_box = hyper_box(
        fill((fix_pts[end] + fix_pts[1]) / 2, D),
        Matrix(Diagonal(fill((fix_pts[end] - fix_pts[1]) / 2, D))),
    )

    if isnothing(name)
        name = "$(space_dim)d-cantor-dust"
    end

    return SelfAffineSet(ifs, fill(T(1//L^D), L^D), b_ball, b_box, name)
end
