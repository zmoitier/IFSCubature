"""Affine transformation."""
struct AffineTransformation{D,T,N}
    A::SMatrix{D,D,T,N}
    b::SVector{D,T}
end

"""Return the affine transformation of g ∘ f."""
function ∘(
    g::AffineTransformation{D,T,N}, f::AffineTransformation{D,T,N}
)::AffineTransformation{D,T,N} where {D,T<:Real,N}
    return AffineTransformation(g.A * f.A, g.A * f.b + g.b)
end

"""Return the translation."""
function translation(vector::Vector{T})::AffineTransformation where {T<:Real}
    D = length(vector)
    N = D * D

    return AffineTransformation{D,T,N}(SMatrix{D,D,T,N}(I), SVector{D,T}(vector))
end

"""Homothety of ratio = `ratio` in dimension `dim`."""
function homothety(ratio::T, D::Int)::AffineTransformation where {T<:Real}
    N = D * D

    return AffineTransformation{D,T,N}(SMatrix{D,D,T,N}(ratio * I), SVector{D,T}(zeros(D)))
end

"""Reflection with respect to the hyperplane define by its normal."""
function reflection(normal::Vector{T})::AffineTransformation where {T<:Real}
    D = length(normal)
    N = D * D

    n = SVector{D,T}(normal)
    t = transpose(n)

    norm2 = t * n
    if isapprox(norm2, 0)
        throw(DomainError(normal, "must be a non-zero vector."))
    end

    return AffineTransformation{D,T,N}(
        SMatrix{D,D,T,N}(I) - (2 / norm2) * (n * t), SVector{D,T}(zeros(D))
    )
end

"""2d Rotation of angle = `angle`."""
function rotation_2d(
    angle::T; implicit_pi::Bool=false
)::AffineTransformation{2,T,4} where {T<:Real}
    A = SMatrix{2,2,T,4}(rotation_matrix_2d(angle; implicit_pi=implicit_pi))
    return AffineTransformation(A, SVector{2,T}(zeros(2)))
end

"""Return the 2d rotation matrix of angle `angle`."""
function rotation_matrix_2d(angle::T; implicit_pi::Bool=false)::Matrix{T} where {T<:Real}
    if implicit_pi
        s, c = sincospi(angle)
    else
        s, c = sincos(angle)
    end

    return [c -s; s c]
end

"""3d Rotation of angle = `angle`."""
function rotation_3d(
    angle::T, axis::Vector{T}; implicit_pi::Bool=false
)::AffineTransformation{3,T,9} where {T<:Real}
    A = SMatrix{3,3,T,9}(rotation_matrix_3d(angle, axis; implicit_pi=implicit_pi))
    return AffineTransformation(A, SVector{3,T}(zeros(3)))
end

"""Return the 3d rotation matrix of angle `angle`."""
function rotation_matrix_3d(
    angle::T, axis::Vector{T}; implicit_pi::Bool=false
)::Matrix{T} where {T<:Real}
    n = norm(axis)

    if isapprox(n, 0)
        throw(DomainError(axis, "must be a non-zero vector."))
    end

    if implicit_pi
        s, c = sincospi(angle)
    else
        s, c = sincos(angle)
    end

    u = axis / n

    cross_product_matrix = SMatrix{3,3,T,9}([
        0 -u[3] u[2]
        u[3] 0 -u[1]
        -u[2] u[1] 0
    ])

    outer_product = u * u'

    return c * I + s * cross_product_matrix + (1 - c) * outer_product
end

"""Return the affine transformation from the points map x → y."""
function affine_transform_from_x_to_y(
    x::SVector{K,SVector{D,T}}, y::SVector{K,SVector{D,T}}
) where {D,T<:Real,K}
    @assert D + 1 == K "must have D+1=$(D+1) equal K=$K."

    X = ones(T, (D + 1, K))
    X[1:D, :] = reinterpret(reshape, T, x)

    Y = reinterpret(reshape, T, y)

    S = Y / X

    N = D * D
    return AffineTransformation{D,T,N}(SMatrix{D,D,T,N}(S[:, 1:D]), SVector{D,T}(S[:, K]))
end
