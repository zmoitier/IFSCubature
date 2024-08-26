### A Pluto.jl notebook ###
# v0.19.46

using Markdown
using InteractiveUtils

# ╔═╡ fe657a18-63c0-11ef-0959-15170cb4d30a
begin
    using Pkg: Pkg
    Pkg.activate(Base.current_project())
    Pkg.instantiate()

    using Printf, LinearAlgebra, StaticArrays, DataStructures, TOML, CairoMakie

    import IFSCubature as src
end

# ╔═╡ b5d45ad9-3996-401f-96e5-78aba17e8686
begin
    #! Type of points
    # const type_points = "equispaced"
    const type_points = "Chebyshev-1"
    # const type_points = "Chebyshev-2"
    # const type_points = "Gauss-Legendre"
    # const type_points = "Gauss-Lobatto"

    const MAXITER = 2000

    #! Data generation
    const save_data = false # take some time

    #! Ploting constants
    const add_title = false
    const save_plot = true

    "Global parameters"
end

# ╔═╡ 252907f1-b344-4466-8f1a-8ff7ea94d7f7
function green_kernel(k, x, y)
    nxy = norm(x - y)
    return exp(im * k * nxy) / nxy
end

# ╔═╡ f40a28d9-7577-4632-94b7-a617d91c1184
function reference_h(
    attractor::src.Attractor{D,T,N}, fct::Function, k::Int, p::Int
) where {D,T<:Real,N}
    quad = src.compute_quadrature(attractor, "Chebyshev-1", k)
    for _ in 1:(p - 1)
        quad = src.refine(attractor, quad)
    end
    res0 = quad(fct)
    quad = src.refine(attractor, quad)
    res1 = quad(fct)

    rel_err = abs((res0 / res1 - 1))

    return (res1, rel_err)
end

# ╔═╡ c9299045-5429-43ee-a3d1-4a1817f09512
function sequence_h_version(
    attractor::src.Attractor{D,T,N},
    quad::src.Quadrature{D,T},
    fct::Function,
    nb_pts_max::Int,
) where {D,T<:Real,N}
    nb_pts = [length(quad)]
    values = [quad(fct)]

    while length(quad) ≤ nb_pts_max
        quad = src.refine(attractor, quad)

        push!(nb_pts, length(quad))
        push!(values, quad(fct))
    end

    return (nb_pts, values)
end

# ╔═╡ 131d23e1-6451-419f-afce-edd809c7fc67
function sequence_p_version(
    attractor::src.Attractor{D,T,N}, type_points::String, fct::Function, nb_pts_max::Int
) where {D,T<:Real,N}
    p = 2
    quad = src.compute_quadrature(attractor, type_points, p; maxiter=MAXITER)

    nb_pts = [length(quad)]
    values = [quad(fct)]
    p += 1

    while length(quad) ≤ nb_pts_max
        quad = src.compute_quadrature(attractor, type_points, p; maxiter=MAXITER)

        push!(nb_pts, length(quad))
        push!(values, quad(fct))

        p += 1
    end

    return (nb_pts, values)
end

# ╔═╡ 9bbf2b35-a4b1-4943-b519-5149a2be2f1f
function result_to_dict(result)
    n = result[1]
    v = result[2]

    return OrderedDict("nb_pts" => n, "values-real" => real.(v), "values-imag" => imag.(v))
end

# ╔═╡ b8fdbae2-445d-444f-8561-3fc953a84983
function _save_data(;
    attractor::src.Attractor{D,T,N}, k::Real, x0::SVector{D,T}, Nh::Int, Np::Int
) where {D,T<:Real,N}
    fct = x -> green_kernel(k, x, x0)
    point_type = "Chebyshev-1"

    data = OrderedDict("k" => k, "x0" => x0, "point-type" => "Chebyshev-1")

    result, precision = reference_h(attractor, fct, 7, 5)
    data["reference"] = OrderedDict(
        "value-real" => real(result), "value-imag" => imag(result), "precision" => precision
    )

    data["barycenter-rule"] = result_to_dict(
        sequence_h_version(attractor, src.barycenter_rule(attractor), fct, Nh)
    )
    for (i, p) in enumerate(2:6)
        data["h-version-Q$i"] = result_to_dict(
            sequence_h_version(
                attractor, src.compute_quadrature(attractor, point_type, p), fct, Nh
            ),
        )
    end
    data["p-version"] = result_to_dict(sequence_p_version(attractor, point_type, fct, Np))

    open("../data-convergences/$(attractor.name).toml", "w") do io
        TOML.print(io, data)
    end

    return nothing
end

# ╔═╡ ef5ef526-11f9-4a6b-bce1-d21870fea3c1
if save_data
    for attractor in [
        src.vicsek_2d(1 / 3),
        src.vicsek_2d(1 / 3, true),
        src.sierpinski_triangle(0.5),
        src.koch_snowflake(),
    ]
        _save_data(; attractor=attractor, k=5.0, x0=SVector{2}([2.0, 0.5]), Nh=2048, Np=626)
    end
end

# ╔═╡ e7504b5e-edd4-49c2-8227-eaf9019bb8ee
function compute_relative_error(seq::Vector{<:Number}, result::Number)
    ε = eps(real(eltype(seq))) / 10

    if isapprox(result, 0)
        return abs.(seq) .+ ε
    end

    return abs.(seq ./ result .- 1) .+ ε
end

# ╔═╡ 5d5bb600-9b89-48b3-90bc-4da6e9f48920
function vicsek_pv()
    fig = Figure()

    ax_args = Dict(
        :xscale => log10, :xlabel => L"M", :yscale => log10, :ylabel => L"Relative error$$"
    )
    if add_title
        ax_args[:title] = L"$p$-version convergence for the Vicsek"
    end
    ax = Axis(fig[1, 1]; ax_args...)

    x_min, x_max = typemax(Int), 0
    for (s, m, l) in [("", :cross, "without"), ("-rot", :xcross, "with")]
        data = TOML.parsefile("../data-convergences/2d-vicsek$s.toml")

        val_ref = data["reference"]["value-real"] + im .* data["reference"]["value-imag"]

        nb_pts = data["p-version"]["nb_pts"]
        val = data["p-version"]["values-real"] .+ im .* data["p-version"]["values-imag"]

        scatterlines!(
            ax,
            nb_pts,
            compute_relative_error(val, val_ref);
            marker=m,
            linestyle=:dash,
            label=L"%$l rotation$$",
        )

        x_min = min(x_min, nb_pts[1])
        x_max = max(x_max, nb_pts[end])
    end

    limits!(ax, (x_min / 1.1, x_max * 1.1), (1e-16, 10))
    axislegend(ax; position=:rt)

    if save_plot
        save("2d-vicsek-pv.pdf", fig)
    end

    return fig
end

# ╔═╡ 71c6ed36-30c4-45d9-8cbb-d0d3c7212f62
vicsek_pv()

# ╔═╡ 15ecd841-9e9e-4291-900c-da1f6886a005
function vicsek_hv()
    fig = Figure()

    ax_args = Dict(
        :xscale => log10, :xlabel => L"M", :yscale => log10, :ylabel => L"Relative error$$"
    )
    if add_title
        ax_args[:title] = L"$h$-version convergence for the Vicsek"
    end
    ax = Axis(fig[1, 1]; ax_args...)

    x_min, x_max = typemax(Int), 0
    for (s, mk, ls, leg) in [("", :cross, :dash, ""), ("-rot", :xcross, :dot, "rot")]
        data = TOML.parsefile("../data-convergences/2d-vicsek$s.toml")

        val_ref = data["reference"]["value-real"] + im .* data["reference"]["value-imag"]

        for k in 1:5
            name = "h-version-Q$k"
            nb_pts = data[name]["nb_pts"]
            val = data[name]["values-real"] .+ im .* data[name]["values-imag"]

            scatterlines!(
                ax,
                nb_pts,
                compute_relative_error(val, val_ref);
                color=k,
                colormap=:tab10,
                colorrange=(1, 10),
                marker=mk,
                linestyle=ls,
                label=L"$Q%$k$ %$leg",
            )

            x_min = min(x_min, nb_pts[1])
            x_max = max(x_max, nb_pts[end])
        end
    end

    limits!(ax, (x_min / 1.1, x_max * 1.1), (1e-16, 10))
    axislegend(ax; position=:lb)

    if save_plot
        save("2d-vicsek-hv.pdf", fig)
    end

    return fig
end

# ╔═╡ 9356d11c-bf4c-4960-90e6-b787ed55a3e4
vicsek_hv()

# ╔═╡ a0b5d444-1dec-4d63-ab09-bd3ee8e02bc3
function sierpinski_koch_pv()
    fig = Figure()

    ax_args = Dict(
        :xscale => log10, :xlabel => L"M", :yscale => log10, :ylabel => L"Relative error$$"
    )
    if add_title
        ax_args[:title] = L"$p$-version convergence for the Vicsek"
    end
    ax = Axis(fig[1, 1]; ax_args...)

    x_min, x_max = typemax(Int), 0
    for (name, mk, leg) in [
        ("sierpinski-triangle", :cross, "Sierpinski triangle"),
        ("koch-snowflake", :xcross, "Koch snowflake"),
    ]
        data = TOML.parsefile("../data-convergences/2d-$name.toml")

        val_ref = data["reference"]["value-real"] + im .* data["reference"]["value-imag"]

        nb_pts = data["p-version"]["nb_pts"]
        val = data["p-version"]["values-real"] .+ im .* data["p-version"]["values-imag"]

        scatterlines!(
            ax,
            nb_pts,
            compute_relative_error(val, val_ref);
            marker=mk,
            linestyle=:dash,
            label=L"%$leg$$",
        )

        x_min = min(x_min, nb_pts[1])
        x_max = max(x_max, nb_pts[end])
    end

    limits!(ax, (x_min / 1.1, x_max * 1.1), (1e-16, 10))
    axislegend(ax; position=:rt)

    if save_plot
        save("2d-sierpinski-koch-pv.pdf", fig)
    end

    return fig
end

# ╔═╡ 7bd6dcf2-71b8-4cf4-b043-3c101a327e08
sierpinski_koch_pv()

# ╔═╡ 12535b8a-b083-4704-b889-db0a460120ff
function sierpinski_koch_hv()
    fig = Figure()

    ax_args = Dict(
        :xscale => log10, :xlabel => L"M", :yscale => log10, :ylabel => L"Relative error$$"
    )
    if add_title
        ax_args[:title] = L"$p$-version convergence for the Vicsek"
    end
    ax = Axis(fig[1, 1]; ax_args...)

    x_min, x_max = typemax(Int), 0
    for (name, mk, ls, leg) in [
        ("sierpinski-triangle", :cross, :dash, "Sierpinski triangle"),
        ("koch-snowflake", :xcross, :dot, "Koch snowflake"),
    ]
        data = TOML.parsefile("../data-convergences/2d-$name.toml")

        val_ref = data["reference"]["value-real"] + im .* data["reference"]["value-imag"]

        for k in 1:5
            name = "h-version-Q$k"
            nb_pts = data[name]["nb_pts"]
            val = data[name]["values-real"] .+ im .* data[name]["values-imag"]

            scatterlines!(
                ax,
                nb_pts,
                compute_relative_error(val, val_ref);
                color=k,
                colormap=:tab10,
                colorrange=(1, 10),
                marker=mk,
                linestyle=ls,
                label=L"$Q%$k$ %$leg",
            )

            x_min = min(x_min, nb_pts[1])
            x_max = max(x_max, nb_pts[end])
        end
    end

    limits!(ax, (x_min / 1.1, x_max * 1.1), (1e-16, 10))
    axislegend(ax; position=:lb)

    if save_plot
        save("2d-sierpinski-koch-hv.pdf", fig)
    end

    return fig
end

# ╔═╡ 50c9daf2-7ebc-434f-a582-06c1cf6d0e4f
sierpinski_koch_hv()

# ╔═╡ Cell order:
# ╠═fe657a18-63c0-11ef-0959-15170cb4d30a
# ╠═b5d45ad9-3996-401f-96e5-78aba17e8686
# ╠═ef5ef526-11f9-4a6b-bce1-d21870fea3c1
# ╠═5d5bb600-9b89-48b3-90bc-4da6e9f48920
# ╠═71c6ed36-30c4-45d9-8cbb-d0d3c7212f62
# ╠═15ecd841-9e9e-4291-900c-da1f6886a005
# ╠═9356d11c-bf4c-4960-90e6-b787ed55a3e4
# ╠═a0b5d444-1dec-4d63-ab09-bd3ee8e02bc3
# ╠═7bd6dcf2-71b8-4cf4-b043-3c101a327e08
# ╠═12535b8a-b083-4704-b889-db0a460120ff
# ╠═50c9daf2-7ebc-434f-a582-06c1cf6d0e4f
# ╠═b8fdbae2-445d-444f-8561-3fc953a84983
# ╠═252907f1-b344-4466-8f1a-8ff7ea94d7f7
# ╠═f40a28d9-7577-4632-94b7-a617d91c1184
# ╠═c9299045-5429-43ee-a3d1-4a1817f09512
# ╠═131d23e1-6451-419f-afce-edd809c7fc67
# ╠═9bbf2b35-a4b1-4943-b519-5149a2be2f1f
# ╠═e7504b5e-edd4-49c2-8227-eaf9019bb8ee
