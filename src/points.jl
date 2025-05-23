function equispaced_points(nb_points::Int; kind::Int=1)
    if kind == 1 # no end points
        @assert nb_points ≥ 1 "nb_points = $nb_points must be ≥ 1."
        n = nb_points - 1
        return Float64.(range(-n, n; step=2) .// nb_points)
    end

    if kind == 2 # with end points
        @assert nb_points ≥ 2 "nb_points = $nb_points must be ≥ 2."
        n = nb_points - 1
        return Float64.(range(-n, n; step=2) .// (nb_points - 1))
    end

    @error "kind must be 1 or 2."
end

function chebyshev_points(nb_points::Int; kind::Int=1)
    if kind == 1 # no end points
        @assert nb_points ≥ 1 "nb_points = $nb_points must be ≥ 1."
        return cospi.((range(nb_points, 1; step=-1) .- 1//2) .// nb_points)
    end

    if kind == 2 # with end points
        @assert nb_points ≥ 2 "nb_points = $nb_points must be ≥ 2."
        return cospi.((range(nb_points, 1; step=-1) .- 1) .// (nb_points - 1))
    end

    @error "kind must be 1 or 2."
end

function gausslegendre_points(nb_points::Int)
    @assert nb_points ≥ 1 "nb_points = $nb_points must be ≥ 1."

    x, _ = gausslegendre(nb_points)
    return x
end

function gausslobatto_points(nb_points::Int)
    @assert nb_points ≥ 2 "nb_points = $nb_points must be ≥ 2."

    x, _ = gausslobatto(nb_points)
    return x
end
