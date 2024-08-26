#!#############################
#! Helper functions
#!#############################

"""Compute the Haussdorff dimension from the contraction factors."""
function _h_dimension(vec_ρ::Vector{Float64})::Float64
    if length(vec_ρ) == 1
        return 0.0
    end

    bracket = extrema((-log(length(vec_ρ))) ./ log.(vec_ρ))

    if isapprox(bracket[1], bracket[2])
        return bracket[1]
    end

    return find_zero(d -> sum(vec_ρ .^ d) .- 1, bracket, Brent())
end

"""Return the bounding ball."""
function _b_ball(ifs::Vector{ContractiveSimilarity{D,T,N}})::HyperBall{D,T} where {D,T,N}
    L = length(ifs)
    A = sum(S.A for S in ifs) ./ L
    b = sum(S.b for S in ifs) ./ L

    z = (I - A) \ b
    r = maximum([norm(S.A * z + S.b - z) / (1 - S.factor) for S in ifs])

    return HyperBall{D,T}(z, r)
end

"""Warn if the open set condition is not satisfy."""
function _cantor_open_set_condition(
    factors::Vector{<:Real}, fix_pts::Vector{<:Real}
)::Nothing
    p_min, p_max = fix_pts[1], fix_pts[end]
    if !all(
        ρᵢ * (p_max - cᵢ) + cᵢ <= ρⱼ * (p_min - cⱼ) + cⱼ for
        ((ρᵢ, ρⱼ), (cᵢ, cⱼ)) in zip(partition(factors, 2, 1), partition(fix_pts, 2, 1))
    )
        @warn "open set condition not satisfy."
    end

    return nothing
end

#!#############################
#! Non fractal IFS
#!#############################

"""Return the segment [a, b]."""
function segment(a::Real, b::Real)::Attractor
    @assert a ≤ b "must have a ≤ b."
    @assert typeof(a) == typeof(b) "a and b must have the same type."

    T = typeof(a)

    factor = T(1//2)
    ifs = [ContractiveSimilarity(factor, [a]), ContractiveSimilarity(factor, [b])]

    h_dim = 1.0

    b_ball = HyperBall((a + b) / 2, (b - a) / 2)

    b_box = HyperBox((a + b) / 2, (b - a) / 2)

    return Attractor{1,T,1}("1d-segment", ifs, h_dim, b_ball, b_box)
end

"""Return the square [a, b]x[a, b]."""
function square(a::Real, b::Real)::Attractor
    @assert a ≤ b "must have a ≤ b."
    @assert typeof(a) == typeof(b) "a and b must have the same type."

    T = typeof(a)

    factor = T(1//2)
    ifs::Vector{ContractiveSimilarity{2,T,4}} = []
    for c in product((a, b), (a, b))
        push!(ifs, ContractiveSimilarity(factor, [v for v in c]))
    end

    center = (a + b) / 2
    radius = (b - a) / 2
    b_ball = HyperBall(fill(center, 2), sqrt(2) * radius)
    b_box = HyperBox(fill(center, 2), Matrix(Diagonal(fill(radius, 2))))

    return Attractor{2,T,4}("2d-square", ifs, 2, b_ball, b_box)
end

"""Return the triangle."""
function triangle()::Attractor
    ifs = [ContractiveSimilarity(0.5, Matrix(-I, 2, 2), [0.0, 0.0])]
    for c in [[1.0, 0.0], [-0.5, √3 / 2], [-0.5, -√3 / 2]]
        push!(ifs, ContractiveSimilarity(0.5, c))
    end

    b_ball = HyperBall(fill(0.0, 2), 1.0)

    b_box = HyperBox([0.25, 0], Matrix(Diagonal([0.75, √3 / 2])))

    return Attractor{2,Float64,4}("2d-triangle", ifs, 2, b_ball, b_box)
end

"""Return the cube [a, b]x[a, b]x[a, b]."""
function cube(a::Real, b::Real)::Attractor
    @assert a ≤ b "must have a ≤ b."
    @assert typeof(a) == typeof(b) "a and b must have the same type."

    T = typeof(a)

    return cantor_dust(T(1//2), [a, b], 3)

    T = typeof(a)

    factor = T(1//2)
    ifs::Vector{ContractiveSimilarity{3,T,9}} = []
    for c in product((a, b), (a, b), (a, b))
        push!(ifs, ContractiveSimilarity(factor, [v for v in c]))
    end

    center = (a + b) / 2
    radius = (b - a) / 2
    b_ball = HyperBall(fill(center, 3), sqrt(3) * radius)
    b_box = HyperBox(fill(center, 3), Matrix(Diagonal(fill(radius, 3))))

    return Attractor{3,T,9}("3d-square", ifs, 3, b_ball, b_box)
end

"""Return the tetrahedron."""
function tetrahedron()::Attractor
    Rz = rotation_matrix_3d(0.0, [1.0, 1.0, 1.0]; implicit_pi=true)

    ifs::Vector{ContractiveSimilarity{3,Float64,9}} = [
        ContractiveSimilarity(0.9, -Rz, fill(1 / 4, 3))
    ]
    for c in [[0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1]]
        push!(ifs, ContractiveSimilarity(0.5, c))
    end

    b_ball = _b_ball(ifs)

    b_box = HyperBox(fill(0.5, 3), Matrix(0.5 * I, 3, 3))

    return Attractor{3,Float64,9}("3d-tetrahedron", ifs, 3, b_ball, b_box)
end

#!#############################
#! Fractal IFS
#!#############################

"""Return the Cantor set."""
function cantor_set(
    contraction_factors::Vector{<:Real}, fix_points::Vector{<:Real}
)::Attractor
    T = eltype(contraction_factors)
    L = length(fix_points)
    idx_sort = sortperm(fix_points)

    @assert all(0 .≤ contraction_factors .< 1) "the contraction factors must be between [0, 1)."
    @assert L ≥ 1 "There must be at least 1 fix points."
    @assert L == length(contraction_factors) "contraction_factors and fix_points must have the same length."

    _cantor_open_set_condition(contraction_factors[idx_sort], fix_points[idx_sort])

    ifs = [
        ContractiveSimilarity(ρ, [c]) for
        (ρ, c) in zip(view(contraction_factors, idx_sort), view(fix_points, idx_sort))
    ]

    h_dim = _h_dimension(contraction_factors)

    b_ball = HyperBall(
        (fix_points[idx_sort[end]] + fix_points[idx_sort[1]]) / 2,
        (fix_points[idx_sort[end]] - fix_points[idx_sort[1]]) / 2,
    )

    b_box = HyperBox(
        (fix_points[idx_sort[end]] + fix_points[idx_sort[1]]) / 2,
        (fix_points[idx_sort[end]] - fix_points[idx_sort[1]]) / 2,
    )

    return Attractor{1,T,1}("1d-cantor-set", ifs, h_dim, b_ball, b_box)
end

"""Return the Cantor set."""
function cantor_set(contraction_factor::Real, fix_points::Vector{<:Real})::Attractor
    @assert 0 ≤ contraction_factor ≤ 1 "the contraction factor must be between [0, 1)."

    return cantor_set(fill(contraction_factor, length(fix_points)), fix_points)
end

"""Return the cantor dust for dimension ≥ 2."""
function cantor_dust(
    contraction_factor::Real, fix_points_1d::Vector{<:Real}, space_dim::Int
)::Attractor
    @assert 0 ≤ contraction_factor ≤ 1 "the contraction factor must be between [0, 1)."
    @assert length(fix_points_1d) ≥ 1 "There must be at least 1 fix points."
    @assert space_dim ≥ 2 "the dimension must be greater or equal than 2."

    D = space_dim
    T = typeof(contraction_factor)
    N = D * D

    fix_pts = sort(fix_points_1d)

    _cantor_open_set_condition(fill(contraction_factor, length(fix_points_1d)), fix_pts)

    ifs::Vector{ContractiveSimilarity{D,T,N}} = []
    for c in product([fix_pts for _ in 1:space_dim]...)
        push!(ifs, ContractiveSimilarity(contraction_factor, [v for v in c]))
    end

    h_dim = -D * log(length(fix_pts)) / log(contraction_factor)

    b_ball = HyperBall(
        fill((fix_pts[end] + fix_pts[1]) / 2, D), sqrt(D) * (fix_pts[end] - fix_pts[1]) / 2
    )

    b_box = HyperBox(
        fill((fix_pts[end] + fix_pts[1]) / 2, D),
        Matrix(Diagonal(fill((fix_pts[end] - fix_pts[1]) / 2, D))),
    )

    return Attractor{D,T,N}("$(space_dim)d-cantor-dust", ifs, h_dim, b_ball, b_box)
end

"""Return the Sierpinski triangle."""
function sierpinski_triangle(contraction_factor::Union{Real,Nothing}=nothing)::Attractor
    if isnothing(contraction_factor)
        factor = 0.5
    else
        factor = contraction_factor
    end

    @assert 0 ≤ factor ≤ 1 "the contraction factor must be between [0, 1)."
    if factor > 0.5
        @warn "open set condition not satisfy."
    end

    D = 2
    T = typeof(factor)
    N = D * D

    ifs = [
        ContractiveSimilarity(factor, c) for
        c in [[1.0, 0.0], [-0.5, √3 / 2], [-0.5, -√3 / 2]]
    ]

    h_dim = -log(3) / log(factor)

    b_ball = HyperBall(fill(0.0, 2), 1.0)

    b_box = HyperBox([0.25, 0], Matrix(Diagonal([0.75, √3 / 2])))

    return Attractor{D,T,N}("2d-sierpinski-triangle", ifs, h_dim, b_ball, b_box)
end

"""Return the bad Sierpinski triangle."""
function bad_sierpinski_triangle(
    contraction_factors::NTuple{3,Real}, fix_points::NTuple{3,Vector{<:Real}}
)::Attractor
    D = 2
    T = Float64
    N = D * D

    ifs = [
        ContractiveSimilarity(
            contraction_factors[1], rotation_matrix_2d(2.0), fix_points[1]
        ),
        ContractiveSimilarity(contraction_factors[2], fix_points[2]),
        ContractiveSimilarity(contraction_factors[3], fix_points[3]),
    ]

    h_dim = _h_dimension([contraction_factors...])

    b_ball = _b_ball(ifs)

    b_box = HyperBox([0.25, 0], Matrix(Diagonal([0.75, √3 / 2])))

    return Attractor{D,T,N}("2d-sierpinski-triangle-bad", ifs, h_dim, b_ball, b_box)
end

"""Return the 2d Vicsek fractal."""
function vicsek_2d(
    contraction_factor::Union{Real,Nothing}=nothing, rotation::Bool=false
)::Attractor
    D = 2
    T = typeof(contraction_factor)
    N = D * D
    name = "2d-vicsek"

    r_max = rotation ? 2 / (√T(2) + 4) : T(1//3)

    if isnothing(contraction_factor)
        factor = r_max
    else
        factor = contraction_factor
    end

    @assert 0 ≤ contraction_factor ≤ 1 "the contraction factor must be between [0, 1)."

    if rotation
        Ri = rotation_matrix_2d(0.25; implicit_pi=true)
        name *= "-rot"
    else
        Ri = Matrix(I, 2, 2)
    end

    if contraction_factor > r_max
        @warn "open set condition not satisfy."
    end

    ifs::Vector{ContractiveSimilarity{D,T,N}} = [
        ContractiveSimilarity(contraction_factor, Ri, [0, 0])
    ]
    for c in product([-1, 1], [-1, 1])
        push!(ifs, ContractiveSimilarity(contraction_factor, [v for v in c]))
    end

    h_dim = -log(5) / log(contraction_factor)

    b_ball = _b_ball(ifs)

    b_box = HyperBox(fill(0.0, 2), Matrix(Diagonal(fill(1.0, 2))))

    return Attractor{D,T,N}(name, ifs, h_dim, b_ball, b_box)
end

"""Return the Sierpinski carpet."""
function sierpinski_carpet(contraction_factor::Union{Real,Nothing}=nothing)::Attractor
    if isnothing(contraction_factor)
        factor = 1 / 3
    else
        factor = contraction_factor
    end

    D = 2
    T = typeof(factor)
    N = D * D

    ifs = [
        ContractiveSimilarity(factor, c) for
        c in [[1, 0], [1, 1], [0, 1], [-1, 1], [-1, 0], [-1, -1], [0, -1], [1, -1]]
    ]

    h_dim = log(8) / log(3)

    b_ball = HyperBall(fill(0.0, D), sqrt(2))

    b_box = HyperBox(fill(0.0, 2), Matrix(Diagonal(fill(1.0, 2))))

    return Attractor{D,T,N}("2d-sierpinski-carpet", ifs, h_dim, b_ball, b_box)
end

"""Return the Koch snowflake."""
function koch_snowflake()::Attractor
    R = rotation_matrix_2d(1 / 6; implicit_pi=true)

    ifs = [ContractiveSimilarity(1 / √3, R, [0, 0])]

    factor = 1 / 3
    for (s, c) in sincospi.((0:5) .// 3)
        push!(ifs, ContractiveSimilarity(factor, [c, s]))
    end

    h_dim = 2.0

    b_ball = _b_ball(ifs)

    b_box = HyperBox(fill(0.0, 2), Matrix(Diagonal([1.0, √3 / 2])))

    D = 2
    T = typeof(factor)
    N = D * D
    return Attractor{D,T,N}("2d-koch-snowflake", ifs, h_dim, b_ball, b_box)
end

"""Return the Gosper flowsnake."""
function gosper_flowsnake()::Attractor
    factor = 1 / √7
    R = rotation_matrix_2d(atan(√3 / 5))

    ifs = [ContractiveSimilarity(factor, R, [0, 0])]
    for (s, c) in sincospi.((0:5) .// 3)
        push!(ifs, ContractiveSimilarity(factor, R, [c, s]))
    end

    h_dim = 2.0

    b_ball = _b_ball(ifs)

    b_box = HyperBox(fill(0.0, 2), Matrix(Diagonal(fill(b_ball.radius, 2))))

    D = 2
    T = typeof(factor)
    N = D * D
    return Attractor{D,T,N}("2d-gosper-flowsnake", ifs, h_dim, b_ball, b_box)
end

"""Return the Sierpinski tetrahedron."""
function sierpinski_tetrahedron(contraction_factor::Union{Real,Nothing}=nothing)::Attractor
    if isnothing(contraction_factor)
        factor = 0.5
    else
        factor = contraction_factor
    end

    @assert 0 ≤ factor ≤ 1 "the contraction factor must be between [0, 1)."

    if factor > 0.5
        @warn "open set condition not satisfy."
    end

    ifs = [
        ContractiveSimilarity(factor, c) for
        c in [[1, 0, -√2 / 2], [-1, 0, -√2 / 2], [0, 1, √2 / 2], [0, -1, √2 / 2]]
    ]

    h_dim = -log(4) / log(factor)

    b_ball = HyperBall(fill(0.0, 3), √(3 / 2))

    b_box = HyperBox(fill(0.0, 3), Matrix(Diagonal([1.0, 1.0, √2 / 2])))

    D = 3
    T = typeof(factor)
    N = D * D
    return Attractor{D,T,N}("3d-sierpinski-tetrahedron", ifs, h_dim, b_ball, b_box)
end

"""Return Menger sponge."""
function menger_sponge(contraction_factor::Union{Real,Nothing}=nothing)::Attractor
    if isnothing(contraction_factor)
        factor = 1 / 3
    else
        factor = contraction_factor
    end

    if factor > 1 / 3
        @warn "open set condition not satisfy."
    end

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

    ifs = [ContractiveSimilarity(factor, c) for c in fix_points]

    h_dim = -log(20) / log(factor)

    b_ball = HyperBall(fill(0.0, 3), sqrt(3))

    b_box = HyperBox(fill(0.0, 3), Matrix(Diagonal(fill(1.0, 3))))

    D = 3
    T = typeof(factor)
    N = D * D
    return Attractor{D,T,N}("3d-menger-sponge", ifs, h_dim, b_ball, b_box)
end

"""Return the 3d Vicsek fractal."""
function vicsek_3d(contraction_factor::Real, central_rotation::Bool=false)::Attractor
    @assert 0 ≤ contraction_factor ≤ 1 "the contraction factor must be between [0, 1)."

    D = 3
    T = typeof(contraction_factor)
    N = D * D

    name = "3d-vicsek"
    if central_rotation
        Ry = rotation_matrix_3d(atan(√2), [0.0, 1.0, 0.0])
        Rz = rotation_matrix_3d(-0.25, [0.0, 0.0, 1.0]; implicit_pi=true)
        R = Rz * Ry
        name *= "-rot"
    else
        R = Matrix(I, 3, 3)
    end
    S0 = ContractiveSimilarity(contraction_factor, R, [0, 0, 0])

    ifs::Vector{ContractiveSimilarity{D,T,N}} = [S0]
    for c in product([-1, 1], [-1, 1], [-1, 1])
        push!(ifs, ContractiveSimilarity(contraction_factor, [v for v in c]))
    end

    h_dim = -log(9) / log(contraction_factor)

    b_ball = HyperBall(fill(0.0, 3), sqrt(3))

    b_box = HyperBox(fill(0.0, 3), Matrix(Diagonal(fill(1.0, 3))))

    return Attractor{D,T,N}(name, ifs, h_dim, b_ball, b_box)
end
