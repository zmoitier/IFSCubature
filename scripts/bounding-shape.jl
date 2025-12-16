using IFSCubature:
    AffineMap,
    HyperBall,
    HyperBox,
    SelfAffineSet,
    affine_map,
    contractive_similarity,
    fix_point,
    fix_points,
    hyper_box_from_corners,
    matrix_rotation_2d,
    similarity_dimension
using LinearAlgebra: I, dot, norm, opnorm
using Optim
using Printf: Format, format
using StaticArrays

function combine(
    ifs_a::Vector{AffineMap{D,T,N}}, ifs_b::Vector{AffineMap{D,T,N}}
) where {D,T,N}
    ifs_c = Vector{AffineMap{D,T,N}}()
    for (f, g) in Iterators.product(ifs_a, ifs_b)
        A, b = f.A * g.A, f.A * g.b + f.b
        push!(ifs_c, AffineMap(A, b, opnorm(A)))
    end
    return ifs_c
end

function smallest_radius(
    z::Union{SVector{D,T},MVector{D,T}}, ifs::Vector{AffineMap{D,T,N}}, fs::Vector{T}
) where {D,T,N}
    return maximum(f * norm(S(z) - z) for (S, f) in zip(ifs, fs))
end

function approximate_bounding_ball(ifs::Vector{AffineMap{D,T,N}}, k::Int=1) where {D,T,N}
    @info "approximate bounding ball"

    fmt = Format("[" * join(fill("%.4f", D), ", ") * "], %.4f")

    ams = deepcopy(ifs)
    ins = 1 ./ (1 .- [S.ρ for S in ams])
    z = zero(SVector{D,T})

    r = smallest_radius(z, ams, ins)
    println(format(fmt, z..., r))

    for _ in 2:k
        ams = combine(ifs, ams)
        ins = 1 ./ (1 .- [S.ρ for S in ams])

        r = smallest_radius(z, ams, ins)
        println(format(fmt, z..., r))
    end

    println()
    return nothing
end

function optimize_bounding_ball(ifs::Vector{AffineMap{D,T,N}}, k::Int=1) where {D,T,N}
    @info "optimize bounding ball"

    fmt = Format("[" * join(fill("%.4f", D), ", ") * "], %.4f")

    ams = deepcopy(ifs)
    ins = 1 ./ (1 .- [S.ρ for S in ams])
    z = MVector{D,T}(sum(fix_point.(ifs)) ./ length(ifs))

    result = Optim.optimize(c -> smallest_radius(c, ams, ins), z)
    # display(result)
    z, r = Optim.minimizer(result), Optim.minimum(result)
    println(format(fmt, z..., r))

    for _ in 2:k
        ams = combine(ifs, ams)
        ins = 1 ./ (1 .- [S.ρ for S in ams])

        result = Optim.optimize(c -> smallest_radius(c, ams, ins), z)
        z, r = Optim.minimizer(result), Optim.minimum(result)
        println(format(fmt, z..., r))
    end

    println()
    return nothing
end

function approximate_bounding_box(
    ifs::Vector{AffineMap{D,T,N}}, measure::Vector{T}, pts_chaos_nb::Int
) where {D,T,N}
    @info "approximate bounding box"

    fmt = Format(
        "[" * join(fill("%.4f", D), ", ") * "], [" * join(fill("%.4f", D), ", ") * "]"
    )

    corner_lower = SVector{D,T}(fill(typemax(T), D))
    corner_upper = SVector{D,T}(fill(typemin(T), D))

    x = zero(SVector{D,T})
    p = cumsum(measure)

    n = pts_chaos_nb ÷ (2 * length(ifs))
    for c in fix_point.(ifs)
        x = c
        corner_lower = min.(corner_lower, x)
        corner_upper = max.(corner_upper, x)

        for _ in 1:n
            k = searchsortedfirst(p, rand())
            x = ifs[k](x)
            corner_lower = min.(corner_lower, x)
            corner_upper = max.(corner_upper, x)
        end

        println(format(fmt, corner_lower..., corner_upper...))

        for _ in 1:n
            k = searchsortedfirst(p, rand())
            x = ifs[k](x)
            corner_lower = min.(corner_lower, x)
            corner_upper = max.(corner_upper, x)
        end

        println(format(fmt, corner_lower..., corner_upper...))
    end

    println()
    return nothing
end

function fct_in_ball(ball::HyperBall{D,T}) where {D,T}
    c, r = ball.center, ball.radius

    function fct(x::SVector{D,T})
        return r - norm(x - c)
    end

    return fct
end

function fct_in_box(box::HyperBox{D,T,N}) where {D,T,N}
    c = box.center
    lengths = SVector{D,T}([norm(v) for v in eachcol(box.paxis)])
    basis = SVector{D,SVector{D,T}}([v / n for (v, n) in zip(eachcol(box.paxis), lengths)])

    function fct(x::SVector{D,T})
        y = x - c
        low, upp = zero(MVector{D,T}), zero(MVector{D,T})
        for (i, (n, e)) in enumerate(zip(lengths, basis))
            p = dot(y, e)
            low[i] = p + n
            upp[i] = n - p
        end
        return (low, upp)
    end

    return fct
end

function check_bounding(
    ifs::Vector{AffineMap{D,T,N}},
    ball::HyperBall{D,T},
    box::HyperBox{D,T,N},
    measure::Vector{T},
    pts_chaos_nb::Int,
) where {D,T,N}
    in_ball = fct_in_ball(ball)
    in_box = fct_in_box(box)

    x = zero(SVector{D,T})
    p = cumsum(measure)

    r_min = typemax(T)
    corner_low = MVector{D,T}(fill(typemax(T), D))
    corner_upp = MVector{D,T}(fill(typemax(T), D))
    for fix_pt in fix_point.(ifs)
        x = fix_pt

        r_min = min(r_min, in_ball(x))
        low, upp = in_box(x)
        corner_low = min.(corner_low, low)
        corner_upp = min.(corner_upp, upp)

        for _ in 1:(pts_chaos_nb ÷ length(ifs))
            k = searchsortedfirst(p, rand())
            x = ifs[k](x)

            r_min = min(r_min, in_ball(x))
            low, upp = in_box(x)
            corner_low = min.(corner_low, low)
            corner_upp = min.(corner_upp, upp)
        end
    end

    @info "Check ball and box: ok" r_min corner_low corner_upp

    return nothing
end

function approximate_stable_bounding_box(
    ifs::Vector{AffineMap{D,T,N}}, n::Int
) where {D,T,N}
    @info "approximate stable bounding box"

    fmt = Format(
        "[" * join(fill("%.4f", D), ", ") * "], [" * join(fill("%.4f", D), ", ") * "]"
    )

    _min = SVector{D,T}(fill(typemax(T), D))
    _max = SVector{D,T}(fill(typemin(T), D))

    for c in fix_point.(ifs)
        _min = min.(c, _min)
        _max = max.(c, _max)
    end

    println(format(fmt, _min..., _max...))

    for _ in 1:n
        for S in ifs
            for x in (_min, SVector(_min[1], _max[2]), SVector(_max[1], _min[2]), _max)
                y = S(x)
                _min = min.(y, _min)
                _max = max.(y, _max)
            end
        end
        println(format(fmt, _min..., _max...))
    end

    println()
    return nothing
end

function brick_2d()
    ρ = 1 / 3
    ifs = [
        contractive_similarity(ρ, c) for
        c in [[-1, -1], [1, -1], [0, 0], [1, 0], [2, 0], [-1, 1], [0, 1], [1, 1], [0, 2]]
    ]
    measure = fill(1 / 9, 9)

    optimize_bounding_ball(ifs, 5)
    approximate_bounding_box(ifs, measure, 500_000)
    approximate_stable_bounding_box(ifs, 2)

    ball = HyperBall(SVector(0.25, 0.25), 1.77)
    box = hyper_box_from_corners([-1.0, -1], [2.0, 2.0])

    check_bounding(ifs, ball, box, measure, 500_000)

    return SelfAffineSet(ifs, measure, ball, box, "2d-brick")
end

function fudgeflake()
    ρ = 1 / √3
    T = matrix_rotation_2d(1 / 6; implicit_pi=true)

    ifs = [
        contractive_similarity(ρ, T, c) for
        c in [[1, 0], [-1 / 2, √3 / 2], [-1 / 2, -√3 / 2]]
    ]
    measure = fill(1 / 3, 3)

    approximate_bounding_ball(ifs, 10)
    approximate_bounding_box(ifs, measure, 500_000)
    approximate_stable_bounding_box(ifs, 20)

    ball = HyperBall(SVector(0.0, 0.0), 1.24)
    box = hyper_box_from_corners([-1.07, -0.98], [1.19, 1.20])

    check_bounding(ifs, ball, box, measure, 500_000)

    return return SelfAffineSet(ifs, measure, ball, box, "2d-fudgeflake")
end

function heighway_dragon()
    ρ = 1 / √2
    ifs = [
        contractive_similarity(ρ, matrix_rotation_2d(1 / 4; implicit_pi=true), [0, 0]),
        contractive_similarity(
            ρ, matrix_rotation_2d(3 / 4; implicit_pi=true), [3 / 5, 1 / 5]
        ),
    ]
    measure = fill(1 / 2, 2)

    optimize_bounding_ball(ifs, 10)
    approximate_bounding_box(ifs, measure, 500_000)
    approximate_stable_bounding_box(ifs, 10)

    ball = HyperBall(SVector(0.41, 0.09), 0.80)
    box = hyper_box_from_corners([-0.34, -0.34], [1.17, 0.67])

    check_bounding(ifs, ball, box, measure, 500_000)

    return SelfAffineSet(ifs, measure, ball, box, "2d-heighway-dragon")
end

function levy_dragon()
    ρ = 1 / √2
    ifs = [
        contractive_similarity(ρ, matrix_rotation_2d(1 / 4; implicit_pi=true), [-1, 0]),
        contractive_similarity(ρ, matrix_rotation_2d(-1 / 4; implicit_pi=true), [1, 0]),
    ]
    measure = fill(1 / 2, 2)

    optimize_bounding_ball(ifs, 10)
    approximate_bounding_box(ifs, measure, 500_000)
    approximate_stable_bounding_box(ifs, 10)

    ball = HyperBall(SVector(0.0, 0.5), 2.07)
    box = hyper_box_from_corners([-2.0, -0.5], [2.0, 2.0])

    check_bounding(ifs, ball, box, measure, 500_000)

    return SelfAffineSet(ifs, measure, ball, box, "2d-levy-dragon")
end

function terdragon()
    ρ = 1 / √3
    Rs = [matrix_rotation_2d(r; implicit_pi=true) for r in [1 / 6, -1 / 2, 1 / 6]]
    xs = [-1, 0, 1]

    ifs = [contractive_similarity(ρ, R, [x, 0]) for (R, x) in zip(Rs, xs)]
    measure = fill(1 / 3, 3)

    approximate_bounding_ball(ifs, 8)
    approximate_bounding_box(ifs, measure, 500_000)
    approximate_stable_bounding_box(ifs, 10)

    ball = HyperBall(SVector(0.0, 0.0), 1.15)
    box = hyper_box_from_corners([-1.13, -0.65], [1.13, 0.65])

    check_bounding(ifs, ball, box, measure, 500_000)

    return SelfAffineSet(ifs, measure, ball, box, "2d-terdragon")
end

function twindragon()
    ρ = 1 / √2
    R = matrix_rotation_2d(1 / 4; implicit_pi=true)

    ifs = [contractive_similarity(ρ, R, [x, 0]) for x in [-1, 1]]
    measure = fill(1 / 2, 2)

    approximate_bounding_ball(ifs, 10)
    approximate_bounding_box(ifs, measure, 500_000)
    approximate_stable_bounding_box(ifs, 10)

    ball = HyperBall(SVector(0.0, 0.0), 1.80)
    box = hyper_box_from_corners([-1.67, -1.34], [1.67, 1.34])

    check_bounding(ifs, ball, box, measure, 500_000)

    return SelfAffineSet(ifs, measure, ball, box, "2d-twindragon")
end

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

    d = similarity_dimension(ifs)
    measure = [S.ρ for S in ifs] .^ d

    optimize_bounding_ball(ifs, 5)
    approximate_bounding_box(ifs, measure, 500_000)
    approximate_stable_bounding_box(ifs, 10)

    ball = HyperBall(SVector(-0.10, 0.10), 1.79)
    box = hyper_box_from_corners([-1.58, -1.11], [1.28, 1.31])

    check_bounding(ifs, ball, box, measure, 500_000)

    return SelfAffineSet(ifs, measure, ball, box, "2d-cantor-non-sym")
end

function barnsley_fern()
    ifs = [
        affine_map([0.0 0.0; 0.0 0.16], [0.0, 0.0]),
        affine_map([0.85 0.04; -0.04 0.85], [0.0, 1.6]),
        affine_map([0.2 -0.26; 0.23 0.22], [0.0, 1.6]),
        affine_map([-0.15 0.28; 0.26 0.24], [0.0, 0.44]),
    ]
    measure = [0.01, 0.85, 0.07, 0.07]

    optimize_bounding_ball(ifs, 5)
    approximate_bounding_box(ifs, measure, 500_000)
    approximate_stable_bounding_box(ifs, 40)

    ball = HyperBall(SVector(1.35, 5.07), 5.25)
    box = hyper_box_from_corners([-2.19, 0.0], [2.66, 10.0])

    check_bounding(ifs, ball, box, measure, 500_000)

    return SelfAffineSet(ifs, measure, ball, box, "2d-barnsley-fern")
end

function brick_3d()
    ρ = 1 / 3
    ifs = [
        contractive_similarity(ρ, c) for c in [
            [-1, -1, -1],
            [0, -1, -1],
            [1, -1, -1],
            [1, 0, -1],
            [1, 1, -1],
            [0, 1, -1],
            [-1, 1, -1],
            [-1, 0, -1],
            #
            [-1, -1, 0],
            [1, -1, 0],
            [1, 0, 0],
            [1, 1, 0],
            [0, 1, 0],
            [-1, 1, 0],
            [0, 0, 0],
            [2, 0, 0],
            [0, 2, 0],
            #
            [-1, -1, 1],
            [0, -1, 1],
            [1, -1, 1],
            [1, 0, 1],
            [1, 1, 1],
            [0, 1, 1],
            [-1, 1, 1],
            [-1, 0, 1],
            [0, 0, 1],
            #
            [0, 0, 2],
        ]
    ]
    measure = fill(1 / 27, 27)

    optimize_bounding_ball(ifs, 3)
    approximate_bounding_box(ifs, measure, 500_000)
    approximate_stable_bounding_box(ifs, 10)

    ball = HyperBall(SVector(0.1, 0.1, 0.1), 1.91)
    box = hyper_box_from_corners([-1.0, -1.0, -1.0], [2.0, 2.0, 2.0])

    check_bounding(ifs, ball, box, measure, 500_000)

    return SelfAffineSet(ifs, measure, ball, box, "3d-brick")
end
