using IFSCubature:
    AffineMap,
    HyperBall,
    HyperBox,
    SelfAffineSet,
    contractive_similarity,
    fix_point,
    fix_points,
    hyper_box_from_corners,
    matrix_rotation_2d
using LinearAlgebra: dot, norm, opnorm, I
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
    c, r2 = ball.center, ball.radius^2

    function fct(x::SVector{D,T})
        return r2 - sum(abs2, x - c)
    end

    return fct
end

function fct_in_box(box::HyperBox{D,T,N}) where {D,T,N}
    c = box.center
    lengths = SVector{D,T}([norm(v) for v in eachcol(box.paxis)])
    basis = SVector{D,SVector{D,T}}([v / n for (v, n) in zip(eachcol(box.paxis), lengths)])

    function fct(x::SVector{D,T})
        y = x - c
        return (
            SVector{D,T}(dot(y, e) + n for (n, e) in zip(lengths, basis)),
            SVector{D,T}(n - dot(y, e) for (n, e) in zip(lengths, basis)),
        )
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

    r = typemax(T)

    for fix_pt in fix_point.(ifs)
        x = fix_pt

        r = min(r, in_ball(x))
        if any(<(0), in_box(x))
            println(in_box(x))
        end

        for _ in 1:(pts_chaos_nb ÷ length(ifs))
            k = searchsortedfirst(p, rand())
            x = ifs[k](x)

            r = min(r, in_ball(x))
            if any(<(0), in_box(x))
                println(in_box(x))
            end
        end
    end

    @info "Check ball and box: ok"
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

    approximate_bounding_ball(ifs, 10)
    approximate_bounding_box(ifs, measure, 500_000)

    ball = HyperBall(SVector(0.0, 0.0), 1.20)
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

    approximate_bounding_ball(ifs, 10)
    approximate_bounding_box(ifs, measure, 500_000)

    ball = HyperBall(SVector(0.0, 0.0), 1.20)
    box = hyper_box_from_corners([-0.34, -0.34], [1.17, 0.67])

    check_bounding(ifs, ball, box, measure, 500_000)

    return SelfAffineSet(ifs, measure, ball, box, "2d-levy-dragon")
end
