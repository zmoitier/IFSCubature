struct LagrangePolynomials{T}
    points::Vector{T}
    weights::Vector{T}
end

function LagrangePolynomials(points::Vector{<:Real})
    N = length(points)
    weights = zero(points)
    for i in 1:N
        weights[i] = 1 / prod(points[i] - points[j] for j in 1:N if j ≠ i)
    end

    return LagrangePolynomials{eltype(points)}(points, weights)
end

length(L::LagrangePolynomials{T}) where {T} = length(L.points)

function (L::LagrangePolynomials{T})(x::Real, i::Int) where {T}
    values = zeros(typeof(x), length(L))
    for (j, point) in enumerate(L.points)
        if isapprox(x, point)
            return typeof(x)(i == j)
        end

        values[j] = (x - point)
    end
    values = L.weights ./ values
    return values[i] / sum(values)
end
