struct Cubature{D,T}
    points::Vector{SVector{D,T}}
    weights::Vector{T}

    function Cubature(points::Vector{SVector{D,T}}, weights::Vector{T}) where {D,T}
        @assert length(points) == length(weights) "The number of points and weights must be the same."
        return new{D,T}(points, weights)
    end
end

length(cbt::Cubature{D,T}) where {D,T} = length(cbt.weights)

function (cbt::Cubature{D,T})(fct::Function) where {D,T}
    return sum(cbt.weights .* fct.(cbt.points))
end

mutable struct HCubature{D,T,N}
    cbt::Cubature{D,T}
    heap::BinaryHeap{Tuple{T,AffineMap{D,T,N},T}} # [ (ρ_w diam Γ, S_w, μ_w) ]

    function HCubature(cbt::Cubature{D,T}, sas::SelfAffineSet{D,T,N}) where {D,T,N}
        return new{D,T,N}(
            cbt,
            BinaryHeap{Tuple{T,AffineMap{D,T,N},T}}(
                Base.Order.By(e -> -e[1]),
                [(
                    2 * sas.bounding_ball.radius,
                    AffineMap(SMatrix{D,D,T,N}(I), SVector{D,T}(zeros(D)), one(T)),
                    one(T),
                )],
            ),
        )
    end
end

function length(hcbt::HCubature{D,T,N}) where {D,T,N}
    return length(hcbt.cbt) * length(hcbt.heap.valtree)
end

function diameter(hcbt::HCubature{D,T,N}) where {D,T,N}
    return first(hcbt.heap)[1]
end

function (hcbt::HCubature{D,T,N})(fct::Function) where {D,T,N}
    s = 0
    for (_, S, μ) in hcbt.heap.valtree
        s += μ * sum(hcbt.cbt.weights .* fct.(S.(hcbt.cbt.points)))
    end
    return s
end

function refine!(hcbt::HCubature{D,T,N}, sas::SelfAffineSet{D,T,N}) where {D,T,N}
    to_refine = [pop!(hcbt.heap)]
    while !isempty(hcbt.heap) && isapprox(first(hcbt.heap)[1], to_refine[1][1])
        push!(to_refine, pop!(hcbt.heap))
    end

    for (hw, Sw, μw) in to_refine
        for (S, μ) in zip(sas.ifs, sas.measure)
            V = Sw ∘ S
            push!(hcbt.heap, (S.ρ * hw, V, μw * μ))
        end
    end

    return nothing
end

function barycenter_rule(sas::SelfAffineSet{D,T,N}) where {D,T,N}
    A = sum([μ .* S.A for (S, μ) in zip(sas.ifs, sas.measure)])
    b = sum([μ .* S.b for (S, μ) in zip(sas.ifs, sas.measure)])

    x0 = (I - A) \ b

    return Cubature([x0], [T(1)])
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
    return Cubature(points, weights)
end
