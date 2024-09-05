struct AffineMap{D,T,N}
    A::SMatrix{D,D,T,N}
    b::SVector{D,T}
    ρ::T
end

function AffineMap(A::SMatrix{D,D,T,N}, b::SVector{D,T}) where {D,T,N}
    return AffineMap(A, b, opnorm(A))
end

function affine_map(
    A::AbstractMatrix, b::AbstractVector, T::Union{DataType,Nothing}=nothing
)
    D = length(b)
    @assert all(x == D for x in size(A)) "The matrix `A` must be a $Dx$D matrix."

    if isnothing(T)
        T = typeof(zero(eltype(A)) + zero(eltype(b)))
    end

    return AffineMap(SMatrix{D,D,T,D * D}(A), SVector{D,T}(b))
end

function ∘(f::AffineMap{D,T,N}, g::AffineMap{D,T,N}) where {D,T,N}
    return AffineMap(f.A * g.A, f.A * g.b + f.b, f.ρ * g.ρ)
end

function (f::AffineMap{1,T,1})(x::T)::T where {T}
    return f.A[1, 1] * x + f.b[1]
end

function (f::AffineMap{D,T,N})(x::Union{SVector{D,T},MVector{D,T}}) where {D,T,N}
    return f.A * x + f.b
end

function (f::AffineMap{D,T,N})(v::Array{T}) where {D,T,N}
    return f.A * v .+ f.b
end

"""Return the affine map from the points map x → y."""
function affine_transform_from_x_to_y(
    x::SVector{K,SVector{D,T}}, y::SVector{K,SVector{D,T}}
) where {D,T,K}
    @assert D + 1 == K "must have D+1=$(D+1) equal K=$K."

    X = ones(T, (D + 1, K))
    X[1:D, :] = reinterpret(reshape, T, x)

    Y = reinterpret(reshape, T, y)

    S = Y / X

    N = D * D
    return AffineMap(SMatrix{D,D,T,N}(S[:, 1:D]), SVector{D,T}(S[:, K]))
end

"""The 2d rotation matrix of angle `angle`."""
function matrix_rotation_2d(angle::Real; implicit_pi::Bool=false)
    if implicit_pi
        s, c = sincospi(angle)
    else
        s, c = sincos(angle)
    end

    return [c -s; s c]
end

"""The 3d rotation matrix of angle `angle`."""
function matrix_rotation_3d(angle::Real, axis::Vector{<:Real}; implicit_pi::Bool=false)
    n = norm(axis)
    @assert !isapprox(n, 0) "`axis` must be a non-zero vector."

    if implicit_pi
        s, c = sincospi(angle)
    else
        s, c = sincos(angle)
    end

    u = axis / n

    cross_product_matrix = [
        0 -u[3] u[2]
        u[3] 0 -u[1]
        -u[2] u[1] 0
    ]

    outer_product = u * u'

    return c * I + s * cross_product_matrix + (1 - c) * outer_product
end

"""The reflection matrix with respect to the hyperplane defined by its normal."""
function matrix_reflection(normal::Vector{<:Real})
    n = normal
    t = transpose(n)

    norm2 = t * n
    @assert !isapprox(norm2, 0) "`normal` must be a non-zero vector."

    return I - (2 / norm2) * (n * t)
end
