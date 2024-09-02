struct Cubature{D,T}
    points::Vector{SVector{D,T}}
    weights::Vector{T}
    h::T
end

function Cubature(
    points::Vector{T}, weights::Vector{<:Real}, h::Real
) where {T<:Union{Real,Vector{<:Real}}}
    @assert length(points) == length(weights) "The length of points and weights must be the same."
    @assert isapprox(sum(weights), 1) "the sum of the weights must = 1."

    S = eltype(weights)

    if eltype(points) <: Real
        D = 1
        pts = [SVector{1,S}(p) for p in points]
    else
        @assert allequal(length.(points)) "All points must have the same length."
        D = length(points[1])
        pts = [SVector{D,S}(p) for p in points]
    end

    return Cubature{D,S}(pts, weights, h)
end

length(cbt::Cubature{D,T}) where {D,T} = length(cbt.weights)

function integral(fct::Function, cbt::Cubature{D,T}) where {D,T}
    return sum(cbt.weights .* fct.(cbt.points))
end

function integral(fct::Function, cbts::Vector{Cubature{D,T}}) where {D,T}
    s = 0
    for i in length(cbts):-1:1
        s += integral(fct, cbts[i])
    end
    return s
end

function refine!(cbts::Vector{Cubature{D,T}}, sas::SelfAffineSet{D,T,N}) where {D,T,N}
    n = length(cbts)
    for (i, cbt) in enumerate(cbts)
        if !(cbt.h ≈ cbts[1].h)
            n = i
        end
    end

    for _ in 1:n
        cbt = popfirst!(cbts)
        for cbt_ref in [
            Cubature{D,T}(S.(cbt.points), μ .* cbt.weights, S.ρ * cbt.h) for
            (S, μ) in zip(sas.ifs, sas.measure)
        ]
            i = searchsortedlast(cbts, cbt_ref; lt=(a, b) -> a.h < b.h, rev=true)
            insert!(cbts, i + 1, cbt_ref)
        end
    end

    return nothing
end

function barycenter_rule(sas::SelfAffineSet{D,T,N}) where {D,T<:Real,N}
    A = sum([μ .* S.A for (S, μ) in zip(sas.ifs, sas.measure)])
    b = sum([μ .* S.b for (S, μ) in zip(sas.ifs, sas.measure)])

    x0 = (I - A) \ b

    return Cubature{D,T}([x0], [T(1)], sas.bounding_ball.radius)
end

function get_points(type_points::String, nb_points::Int)
    if type_points == "Equispaced-1"
        return equispaced_points(nb_points; kind=1)

    elseif type_points == "Equispaced-2"
        return equispaced_points(nb_points; kind=2)

    elseif type_points == "Chebyshev-1"
        return chebyshev_points(nb_points; kind=1)

    elseif type_points == "Chebyshev-2"
        return chebyshev_points(nb_points; kind=2)

    elseif type_points == "Gauss-Legendre"
        return gausslegendre_points(nb_points)

    elseif type_points == "Gauss-Lobatto"
        return gausslobatto_points(nb_points)

    else
        mgs = "possible type of points are:\n"
        for choice in [
            "equispaced-1",
            "equispaced-2",
            "Chebyshev-1",
            "Chebyshev-2",
            "Gauss-Legendre",
            "Gauss-Lobatto",
        ]
            mgs *= "  > $choice\n"
        end
        throw(ArgumentError(mgs))
    end
end

function matrix_from_lagrange_polynomials(
    sas::SelfAffineSet{1,T,1}, type_points::String, nb_points::Int
) where {T}
    bounding_box = sas.bounding_box

    points =
        only(bounding_box.paxis) .* get_points(type_points, nb_points) .+
        only(bounding_box.center)

    L = LagrangePolynomials(points)
    n = length(L)

    pts = zero(points)
    M = zeros(T, (n, n))
    for (S, μ) in zip(sas.ifs, sas.measure)
        pts = S.A .* points .+ S.b
        for i in 1:n
            M[i, :] += μ .* L.(pts, i)
        end
    end

    return (M, [SVector{1,T}(p) for p in points])
end

function matrix_from_lagrange_polynomials(
    sas::SelfAffineSet{2,T,4}, type_points::String, nb_points::Int
) where {T}
    bounding_box = sas.bounding_box
    v = get_points(type_points, nb_points)

    L = LagrangePolynomials(v)
    n = length(L)
    n2 = n^2

    pts_ref = zeros(T, (2, n2))
    for (i, p) in enumerate(product(v, v))
        pts_ref[:, i] .= p
    end

    pts = bounding_box.paxis * pts_ref .+ bounding_box.center

    pts_tmp = zero(pts_ref)
    M = zeros(T, (n2, n2))
    for (S, μ) in zip(sas.ifs, sas.measure)
        pts_tmp = bounding_box.paxis \ (S.A * pts .+ S.b .- bounding_box.center)
        for (i, (p, q)) in enumerate(product(1:n, 1:n))
            M[i, :] += μ .* L.(view(pts_tmp, 1, :), p) .* L.(view(pts_tmp, 2, :), q)
        end
    end

    return (M, [SVector{2,T}(p) for p in eachcol(pts)])
end

function matrix_from_lagrange_polynomials(
    sas::SelfAffineSet{3,T,9}, type_points::String, nb_points::Int
) where {T}
    bounding_box = sas.bounding_box
    v = get_points(type_points, nb_points)

    L = LagrangePolynomials(v)
    n = length(L)
    n3 = n^3

    pts_ref = zeros(T, (3, n3))
    for (i, p) in enumerate(product(v, v, v))
        pts_ref[:, i] .= p
    end

    pts = bounding_box.paxis * pts_ref .+ bounding_box.center

    pts_tmp = zero(pts_ref)
    M = zeros(T, (n3, n3))
    for (S, μ) in zip(sas.ifs, sas.measure)
        pts_tmp = bounding_box.paxis \ (S.A * pts .+ S.b .- bounding_box.center)
        for (i, (p, q, r)) in enumerate(product(1:n, 1:n, 1:n))
            M[i, :] +=
                μ .* L.(view(pts_tmp, 1, :), p) .* L.(view(pts_tmp, 2, :), q) .*
                L.(view(pts_tmp, 3, :), r)
        end
    end

    return (M, [SVector{3,T}(p) for p in eachcol(pts)])
end

function compute_cubature(
    sas::SelfAffineSet{D,T,N},
    type_points::String,
    nb_points::Int;
    maxiter::Union{Nothing,Integer}=nothing,
) where {D,T,N}
    M, points = matrix_from_lagrange_polynomials(sas, type_points, nb_points)
    J = length(points)

    if isnothing(maxiter)
        maxiter = max(100, J)
    end

    x0 = rand(T, J)
    λ, weights = powm!(M, x0; maxiter=maxiter, tol=1e-14)

    @assert isapprox(λ, 1) "$λ should be equal to 1."

    weights /= sum(weights)

    @assert isapprox(sum(weights), 1) "the sum of the weights must = 1."
    return Cubature{D,T}(points, weights, sas.bounding_ball.radius)
end
