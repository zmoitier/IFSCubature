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

    ball = hyper_ball(fill(0.0, 2), 1.0)
    box = hyper_box([0.25, 0], Matrix(Diagonal([0.75, √3 / 2])))

    return SelfAffineSet(ifs, fill(1 / 4, 4), ball, box, "2d-triangle")
end

"""Return the Sierpinski triangle."""
function sierpinski_triangle(contraction_factor::Union{Float64,Nothing}=nothing)
    if isnothing(contraction_factor)
        ρ = 0.5
    else
        ρ = contraction_factor
    end

    @assert 0 ≤ ρ < 1 "the contraction factor must be between [0, 1)."

    ifs = [
        contractive_similarity(ρ, c) for c in [[1.0, 0.0], [-0.5, √3 / 2], [-0.5, -√3 / 2]]
    ]

    ball = hyper_ball(fill(0.0, 2), 1.0)
    box = hyper_box([0.25, 0], Matrix(Diagonal([0.75, √3 / 2])))

    return SelfAffineSet(ifs, fill(1 / 3, 3), ball, box, "2d-sierpinski-triangle")
end

"""Return the Sierpinski triangle."""
function sierpinski_triangle_fat(n::Int)
    @assert n ≥ 2 "`n` must be greater or equal to 2."

    e = collect(1:n)
    ρ = find_zero(x -> sum(x .^ e) .- 1, (0.5, 2 / 3), Brent())
    # τ = find_zero(x -> x .^ (n + 1) .- x .+ 1 / 3, (0.3, 0.4), Brent())
    # d = log(τ) / log(ρ)

    ifs = [
        contractive_similarity(ρ, c) for c in [[1.0, 0.0], [-0.5, √3 / 2], [-0.5, -√3 / 2]]
    ]

    ball = hyper_ball(fill(0.0, 2), 1.0)
    box = hyper_box([0.25, 0], Matrix(Diagonal([0.75, √3 / 2])))

    return SelfAffineSet(ifs, fill(1 / 3, 3), ball, box, "2d-sierpinski-triangle-fat")
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

    ball = hyper_ball(fill(0.0, 2), √2)
    box = hyper_box(fill(0.0, 2), Matrix(Diagonal(fill(1.0, 2))))

    return SelfAffineSet(ifs, fill(1 / 5, 5), ball, box, name)
end

"""Return the Sierpinski carpet."""
function sierpinski_carpet()
    ρ = 1 / 3
    ifs = [
        contractive_similarity(ρ, c) for
        c in [[1, 0], [1, 1], [0, 1], [-1, 1], [-1, 0], [-1, -1], [0, -1], [1, -1]]
    ]

    ball = hyper_ball(fill(0.0, 2), √2)
    box = hyper_box(fill(0.0, 2), Matrix(Diagonal(fill(1.0, 2))))

    return SelfAffineSet(ifs, fill(1 / 8, 8), ball, box, "2d-sierpinski-carpet")
end

"""Return the Koch snowflake."""
function koch_snowflake()
    R = matrix_rotation_2d(1 / 6; implicit_pi=true)

    ifs = [contractive_similarity(1 / √3, R, fill(0.0, 2))]

    ρ = 1 / 3
    for (s, c) in sincospi.((0:5) .// 3)
        push!(ifs, contractive_similarity(ρ, [c, s]))
    end

    ball = hyper_ball(fill(0.0, 2), 1.0)
    box = hyper_box(fill(0.0, 2), Matrix(Diagonal([1.0, √3 / 2])))

    μ = [1 / 3, fill(1 / 9, 6)...]

    return SelfAffineSet(ifs, μ, ball, box, "2d-koch-snowflake")
end

"""Return the Gosper flowsnake."""
function gosper_flowsnake()
    ρ = 1 / √7
    R = matrix_rotation_2d(atan(√3 / 5))

    ifs = [contractive_similarity(ρ, R, fill(0.0, 2))]
    for (s, c) in sincospi.((0:5) .// 3)
        push!(ifs, contractive_similarity(ρ, R, [c, s]))
    end

    ball = hyper_ball(fill(0.0, 2), √3 / (√7 - 1))
    box = hyper_box(fill(0.0, 2), Matrix(Diagonal(fill(ball.radius, 2))))

    return SelfAffineSet(ifs, fill(1 / 7, 7), ball, box, "2d-gosper-flowsnake")
end

"""Return the brick."""
function brick_2d()
    ρ = 1 / 3
    ifs = [
        contractive_similarity(ρ, c) for
        c in [[-1, -1], [1, -1], [0, 0], [1, 0], [2, 0], [-1, 1], [0, 1], [1, 1], [0, 2]]
    ]

    ball = bounding_ball(ifs)
    box = bounding_box(ifs)

    return SelfAffineSet(ifs, fill(1 / 9, 9), ball, box, "2d-brick")
end

"""Return the Durer's Pentagon."""
function durer_pentagon()
    ρ = (3 - √5) / 2
    ifs = [contractive_similarity(ρ, -Matrix(I, 2, 2), zeros(2))]
    for (s, c) in sincospi.((2 / 5) .* (0:4))
        push!(ifs, contractive_similarity(ρ, [c, s]))
    end

    ball = bounding_ball(ifs)
    box = bounding_box(ifs)

    return SelfAffineSet(ifs, fill(1 / 6, 6), ball, box, "2d-durer-pentagon")
end

"""Return the Fudgeflake."""
function fudgeflake()
    ρ = 1 / √3
    T = matrix_rotation_2d(1 / 6; implicit_pi=true)

    ifs = [
        contractive_similarity(ρ, T, c) for
        c in [[1, 0], [-1 / 2, √3 / 2], [1 / 2, -√3 / 2]]
    ]

    ball = bounding_ball(ifs)
    box = bounding_box(ifs)
    if ball.radius < box.paxis[1, 1]
        box = hyper_box(ball.center, [ball.radius 0; 0 ball.radius])
    end

    return SelfAffineSet(ifs, fill(1 / 3, 3), ball, box, "2d-fudgeflake")
end

"""Return the Heighway Dragon."""
function heighway_dragon()
    ρ = 1 / √2
    ifs = [
        contractive_similarity(ρ, matrix_rotation_2d(1 / 4; implicit_pi=true), [0, 0]),
        contractive_similarity(
            ρ, matrix_rotation_2d(3 / 4; implicit_pi=true), [3 / 5, 1 / 5]
        ),
    ]

    ball = bounding_ball(ifs)
    box = hyper_box(ball.center, Diagonal(fill(ball.radius, 2)))

    return SelfAffineSet(ifs, fill(1 / 2, 2), ball, box, "2d-heighway-dragon")
end

"""Return the Lévy Dragon."""
function levy_dragon()
    ρ = 1 / √2
    ifs = [
        contractive_similarity(ρ, matrix_rotation_2d(1 / 4; implicit_pi=true), [-1, 0]),
        contractive_similarity(ρ, matrix_rotation_2d(-1 / 4; implicit_pi=true), [1, 0]),
    ]

    ball = bounding_ball(ifs)
    box = hyper_box(ball.center, Diagonal(fill(ball.radius, 2)))

    return SelfAffineSet(ifs, fill(1 / 2, 2), ball, box, "2d-levy-dragon")
end

"""Return the Terdragon."""
function terdragon()
    ρ = 1 / √3
    Rs = [matrix_rotation_2d(r; implicit_pi=true) for r in [1 / 6, -1 / 2, 1 / 6]]
    xs = [-1, 0, 1]

    ifs = [contractive_similarity(ρ, R, [x, 0]) for (R, x) in zip(Rs, xs)]

    ball = bounding_ball(ifs)
    box = bounding_box(ifs)
    if ball.radius < box.paxis[1, 1]
        box = hyper_box(ball.center, [ball.radius 0; 0 ball.radius])
    end

    return SelfAffineSet(ifs, fill(1 / 3, 3), ball, box, "2d-terdragon")
end

"""Return the Twindragon."""
function twindragon()
    ρ = 1 / √2
    R = matrix_rotation_2d(1 / 4; implicit_pi=true)

    ifs = [contractive_similarity(ρ, R, [x, 0]) for x in [-1, 1]]

    ball = bounding_ball(ifs)
    box = hyper_box(ball.center, Diagonal(fill(ball.radius, 2)))

    return SelfAffineSet(ifs, fill(1 / 2, 2), ball, box, "2d-twindragon")
end

"""Return a nonsymmetric Cantor dust."""
function cantor_dust_non_sym()
    rs::Vector{Float64} = [0.25, 0.35, 0.3, 0.4]
    Ms::Vector{Matrix{Float64}} = [
        matrix_rotation_2d(0.4),
        matrix_rotation_2d(0.2),
        matrix_rotation_2d(0.3),
        matrix_rotation_2d(0.1),
    ]
    vs::Vector{Vector{Float64}} = [[-1.4, -1.1], [0.8, -0.7], [1.2, 1.3], [-1.3, 0.9]]

    ifs = [contractive_similarity(r, M, v) for (r, M, v) in zip(rs, Ms, vs)]

    ball = bounding_ball(ifs)
    box = bounding_box(ifs)
    if ball.radius < box.paxis[1, 1]
        box = hyper_box(ball.center, [ball.radius 0; 0 ball.radius])
    end

    d = dimension(ifs)
    measure = [S.ρ for S in ifs] .^ d[2]

    return SelfAffineSet(ifs, measure, ball, box, "2d-cantor-non-sym")
end

"""Return the Barnsley fern."""
function barnsley_fern()
    ifs = [
        affine_map([0.0 0.0; 0.0 0.16], [0.0, 0.0]),
        affine_map([0.85 0.04; -0.04 0.85], [0.0, 1.6]),
        affine_map([0.2 -0.26; 0.23 0.22], [0.0, 1.6]),
        affine_map([-0.15 0.28; 0.26 0.24], [0.0, 0.44]),
    ]

    ball = bounding_ball(ifs)
    box = bounding_box(ifs)
    if ball.radius < box.paxis[1, 1]
        box = hyper_box(ball.center, [ball.radius 0; 0 ball.radius])
    end

    d = dimension(ifs)[1]
    measure = [svdvals(S.A)[end] for S in ifs] .^ d

    return SelfAffineSet(ifs, measure, ball, box, "2d-barnsley-fern")
end
