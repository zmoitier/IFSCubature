function equispaced_points(nb_points::Int; kind::Int=1)::Vector{Float64}
    if kind == 1 # no end points
        @assert nb_points ≥ 1 "nb_points = $nb_points must be ≥ 1."
        return collect(range(-1, 1, nb_points + 2))[2:(end - 1)]
    end

    if kind == 2 # with end points
        @assert nb_points ≥ 2 "nb_points = $nb_points must be ≥ 2."
        return collect(range(-1, 1, nb_points))
    end

    @error "kind must be 1 or 2."
end

function chebyshev_points(nb_points::Int; kind::Int=1)::Vector{Float64}
    if kind == 1 # no end points
        @assert nb_points ≥ 1 "nb_points = $nb_points must be ≥ 1."
        return cospi.((collect(range(nb_points, 1; step=-1)) .- 0.5) ./ nb_points)
    end

    if kind == 2 # with end points
        @assert nb_points ≥ 2 "nb_points = $nb_points must be ≥ 2."
        return cospi.((collect(range(nb_points, 1; step=-1)) .- 1) ./ (nb_points - 1))
    end

    @error "kind must be 1 or 2."
end

function gausslegendre_points(nb_points::Int)::Vector{Float64}
    @assert nb_points ≥ 1 "nb_points = $nb_points must be ≥ 1."

    x, _ = gausslegendre(nb_points)
    return x
end

function gausslobatto_points(nb_points::Int)::Vector{Float64}
    @assert nb_points ≥ 2 "nb_points = $nb_points must be ≥ 2."

    x, _ = gausslobatto(nb_points)
    return x
end
