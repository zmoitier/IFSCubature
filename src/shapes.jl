struct HyperBall{D,T}
    center::SVector{D,T}
    radius::T
end

function HyperBall(center::Vector{<:Real}, radius::Real)
    @assert radius > 0 "radius=$radius > 0."

    D = length(center)
    T = typeof(radius)

    return HyperBall{D,T}(SVector{D,T}(center), radius)
end

function HyperBall(center::Real, radius::Real)
    @assert radius ≥ 0 "paxis must be non-negative."

    T = typeof(radius)

    return HyperBall{1,T}(SVector{1,T}(center), radius)
end

struct HyperBox{D,T,N}
    center::SVector{D,T}
    paxis::SMatrix{D,D,T,N}
end

function HyperBox(center::Vector{<:Real}, paxis::Matrix{<:Real})
    D = length(center)
    T = eltype(paxis)
    N = D * D

    @assert size(paxis) == (D, D) "size(paxis) ≠ ($D, $D)"
    @assert !isapprox(det(paxis), 0) "the column of paxis must form a basis."

    return HyperBox{D,T,N}(SVector{D,T}(center), SMatrix{D,D,T,N}(paxis))
end

function HyperBox(center::Real, half_width::Real)
    @assert half_width ≥ 0 "paxis must be non-negative."

    T = typeof(half_width)

    return HyperBox{1,T,1}(SVector{1,T}(center), SMatrix{1,1,T,1}(half_width))
end

struct Segment{T}
    vertices::SVector{2,T}
end

function Segment(a::Real, b::Real)
    v = [a, b]
    T = eltype(v)
    return Segment{T}(SVector{2,T}(v))
end

struct Polygon{T}
    vertices::Vector{SVector{2,T}}
end

function Polygon(vertices::Vector{Vector{T}}) where {T<:Real}
    return Polygon{T}([SVector{2,T}(v) for v in vertices])
end

struct Polyhedron{T,S}
    vertices::Vector{SVector{3,T}}
    triangles::NTuple{3,Vector{S}}
    edges::NTuple{2,Vector{S}}
end

function Polyhedron(
    vertices::Vector{Vector{T}}, triangles::NTuple{3,Vector{S}}, edges::NTuple{2,Vector{S}}
) where {T<:Real,S<:Integer}
    return Polyhedron{T,S}([SVector{3,T}(v) for v in vertices], _1_to_0.(triangles), edges)
end

function _1_to_0(v::Vector{S})::Vector{S} where {S<:Integer}
    return v .- 1
end
