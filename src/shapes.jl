struct HyperBall{D,T}
    center::SVector{D,T}
    radius::T
end

function vertices(ball::HyperBall{1,T}) where {T}
    return [ball.center[1] - ball.radius, ball.center[1] + ball.radius]
end

function hyper_ball(center::AbstractVector{T}, radius::T) where {T}
    @assert radius > 0 "radius=$radius > 0."

    D = length(center)
    return HyperBall{D,T}(SVector{D,T}(center), radius)
end

function hyper_ball(center::T, radius::T) where {T}
    @assert radius ≥ 0 "paxis must be non-negative."

    return HyperBall{1,T}(SVector{1,T}(center), radius)
end

function (f::AffineMap{D,T,N})(hb::HyperBall{D,T}) where {D,T,N}
    return HyperBall{D,T}(f(hb.center), hb.radius * f.ρ)
end

struct HyperBox{D,T,N}
    center::SVector{D,T}
    paxis::SMatrix{D,D,T,N}
end

function vertices(box::HyperBox{D,T,N}) where {D,T,N}
    pts = SVector{D,T}[]
    for p in product([[-1, 1] for _ in 1:D]...)
        push!(pts, box.center + box.paxis * SVector{D,T}(p))
    end

    return pts
end

function vertices(box::HyperBox{1,T,1}) where {T}
    return [box.center[1] - box.paxis[1], box.center[1] + box.paxis[1]]
end

function hyper_box(center::AbstractVector{T}, paxis::AbstractMatrix{T}) where {T}
    D = length(center)
    @assert size(paxis) == (D, D) "size(paxis) ≠ ($D, $D)"

    N = D * D
    return HyperBox{D,T,N}(SVector{D,T}(center), SMatrix{D,D,T,N}(paxis))
end

function hyper_box(center::T, half_width::T) where {T}
    return HyperBox{1,T,1}(SVector{1,T}(center), SMatrix{1,1,T,1}(half_width))
end

function (f::AffineMap{D,T,N})(hb::HyperBox{D,T,N}) where {D,T,N}
    return HyperBox{D,T,N}(f(hb.center), f.A * hb.paxis)
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

function Polygon(vertices::Vector{Vector{T}}) where {T}
    return Polygon{T}([SVector{2,T}(v) for v in vertices])
end

struct Polyhedron{T,S}
    vertices::Vector{SVector{3,T}}
    triangles::NTuple{3,Vector{S}}
    edges::NTuple{2,Vector{S}}
end

function Polyhedron(
    vertices::Vector{Vector{T}}, triangles::NTuple{3,Vector{S}}, edges::NTuple{2,Vector{S}}
) where {T,S}
    return Polyhedron{T,S}(
        [SVector{3,T}(v) for v in vertices], _shift_down.(triangles), edges
    )
end

function _shift_down(v::Vector{S}) where {S}
    return v .- 1
end
