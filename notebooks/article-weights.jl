### A Pluto.jl notebook ###
# v0.19.46

using Markdown
using InteractiveUtils

# ╔═╡ 80b8650e-6390-11ef-2417-eb2b5b153f95
begin
    using Pkg: Pkg
    Pkg.activate(Base.current_project())
    Pkg.instantiate()

    using Printf, StaticArrays, CairoMakie

    import IFSCubature as src
end

# ╔═╡ c6939b67-42a4-45e1-9c61-33e9c1397d4c
begin
    #! Type of points
    # const type_points = "equispaced"
    const type_points = "Chebyshev-1"
    # const type_points = "Chebyshev-2"
    # const type_points = "Gauss-Legendre"
    # const type_points = "Gauss-Lobatto"

    #! Data generation
    const save_data = false # take some time

    #! Ploting constants
    const add_title = true
    const save_plot = false

    "Global parameters"
end

# ╔═╡ f788e3a3-181f-4676-9463-ed5b53e47223
function _save_weights(
    attractor::src.Attractor, M_max::Int, suffix::Union{String,Nothing}=nothing
)
    filename = attractor.name
    if !isnothing(suffix)
        filename = filename * "-$suffix"
    end

    open("$filename.csv", "w") do file
        write(file, "nb-points,sum-abs\n")
        write(file, "1,$(@sprintf("%.16e", 1))\n")

        for M in 2:M_max
            quad = src.compute_quadrature(attractor, type_points, M; maxiter=2000)

            n = @sprintf("%d", length(quad))
            w = @sprintf("%.16e", sum(abs.(quad.weights)))
            write(file, "$n,$w\n")
        end
    end

    return nothing
end

# ╔═╡ f68756e5-353a-4481-8e03-0796b5384f90
if save_data
    for (factor, suffix) in [(1 / 2, "1o2"), (1 / 3, "1o3"), (1 / 4, "1o4")]
        print((factor, suffix))
        _save_weights(src.cantor_set(factor, [0, 1]), 500, suffix)
        println(": done")
    end

    for attractor in [src.vicsek_2d(1 / 3), src.vicsek_2d(1 / 3, true)]
        _save_weights(attractor, 32)
    end
end

# ╔═╡ 00cf0514-7931-4d7e-a523-ae918eac6de3
function cantor_set_weights(M::Int)
    attractor = src.cantor_set(1 / 3, [0, 1])
    quad = src.compute_quadrature(attractor, type_points, M)
    w_lim = extrema(quad.weights)
    W = maximum(abs.(w_lim))

    dx, dy = 0.025, 0.05 * (w_lim[2] - w_lim[1])
    X = (-dx, 1 + dx)
    Y = (w_lim[1] - dy, w_lim[2] + dy)

    pf = [[
        only(attractor.bounding_box.center) - only(attractor.bounding_box.paxis),
        only(attractor.bounding_box.center) + only(attractor.bounding_box.paxis),
    ]]
    for _ in 1:5
        pf = [S.(part) for S in attractor.ifs for part in pf]
    end

    fig = Figure()

    ax_args = Dict(:xlabel => L"x", :ylabel => L"w")
    if add_title
        ax_args[:title] = L"Cantor set with $\rho=1/3$ and $ M = %$M $"
    end
    ax = Axis(fig[1, 1]; ax_args...)

    limits!(ax, X, Y)

    for p in pf
        poly!(
            ax,
            [(p[1], Y[1]), (p[2], Y[1]), (p[2], Y[2]), (p[1], Y[2])];
            color=(:black, 0.3),
        )
    end

    scatter!(
        ax,
        reinterpret(Float64, quad.points),
        quad.weights;
        color=quad.weights,
        colormap=:vik,
        colorrange=(-W, W),
        strokewidth=1,
        strokecolor=:black,
    )

    if save_plot
        save("1d-cantor-set-weights.pdf", fig)
    end

    return fig
end

# ╔═╡ baa42609-1d76-439a-a82e-da277cbbc754
cantor_set_weights(100)

# ╔═╡ ef2380b5-1ac8-43a7-be6c-b50d29608af5
function vicsek_weights(M::Int)
    attractor = src.vicsek_2d(1 / 3)
    quad = src.compute_quadrature(attractor, type_points, M)
    W = maximum(abs.(quad.weights))

    pf = [[
        SVector{2,Float64}([-1, -1]),
        SVector{2,Float64}([1, -1]),
        SVector{2,Float64}([1, 1]),
        SVector{2,Float64}([-1, 1]),
    ]]
    for _ in 1:5
        pf = [S.(part) for S in attractor.ifs for part in pf]
    end

    fig = Figure()

    ax_args = Dict(:xlabel => L"x", :ylabel => L"y", :aspect => 1)
    if add_title
        ax_args[:title] = L"Vicsek with $\rho=1/3$ and $ M = %$(M*M) $"
    end
    ax = Axis(fig[1, 1]; ax_args...)

    dx = 0.025
    limits!(ax, (-1 - dx, 1 + dx), (-1 - dx, 1 + dx))

    for part in pf
        poly!(ax, part; color=(:black, 0.5))
    end

    sc = scatter!(
        ax,
        quad.points;
        color=quad.weights,
        colormap=:vik,
        colorrange=(-W, W),
        strokewidth=1,
        strokecolor=:black,
    )

    Colorbar(fig[1, 2], sc; label=L"w")

    if save_plot
        save("2d-vicsek-weights.pdf", fig)
    end

    return fig
end

# ╔═╡ 021d4f9e-5faa-4c22-9d2f-9ba0916d4b73
vicsek_weights(10)

# ╔═╡ 7dbc6f7f-db53-45d6-8378-4dae562608f5
function cantor_set_sum_weights()
    fig = Figure()

    ax_args = Dict(:xlabel => L"M", :ylabel => L"|\mathbf{w}|_1")
    if add_title
        ax_args[:title] = L"Sum of the absolute value of the weight for Cantor set$$"
    end
    ax = Axis(fig[1, 1]; ax_args...)

    x = collect(1:500)
    y = 1 .+ 0.14 .* log.(x)
    lines!(ax, x, y; color=:black, label=L"1 + \mathcal{O}(\log M)")

    for s in [4, 3, 2]
        nb_pts = Int[]
        sum_abs = Float64[]

        open("../data-weights/1d-cantor-set-1o$s.csv", "r") do file
            _ = readline(file)
            for line in eachline(file)
                numbers = split(line, ",")
                push!(nb_pts, parse(Int, numbers[1]))
                push!(sum_abs, parse(Float64, numbers[2]))
            end
        end

        lines!(ax, nb_pts, sum_abs; label=L"\rho = 1/%$s")
    end

    axislegend(ax; position=:lt)

    if save_plot
        save("1d-cantor-set-sum-weights.pdf", fig)
    end

    return fig
end

# ╔═╡ f89e7cb5-f6f7-4fe0-adc0-96b10963e755
cantor_set_sum_weights()

# ╔═╡ a2e84401-9a34-4446-850b-74e729bc82b9
function vicsek_sum_weights()
    fig = Figure()

    ax_args = Dict(:xlabel => L"M", :ylabel => L"|\mathbf{w}|_1")
    if add_title
        ax_args[:title] = L"Sum of the absolute value of the weight for Vicsek$$"
    end
    ax = Axis(fig[1, 1]; ax_args...)

    x = range(1, 1024; length=64)
    y = 1 .+ 0.045 .* log.(x)
    lines!(ax, x, y; color=:black, label=L"1 + \mathcal{O}(\log M)")

    for (s, m) in [("", :cross), ("-rot", :xcross)]
        nb_pts = Int[]
        sum_abs = Float64[]

        open("../data-weights/2d-vicsek$s.csv", "r") do file
            _ = readline(file)
            for line in eachline(file)
                numbers = split(line, ",")
                push!(nb_pts, parse(Int, numbers[1]))
                push!(sum_abs, parse(Float64, numbers[2]))
            end
        end

        scatterlines!(
            ax, nb_pts, sum_abs; marker=m, linestyle=:dash, label=L"without rotation$$"
        )
    end

    axislegend(ax; position=:lt)

    if save_plot
        save("2d-vicsek-sum-weights.pdf", fig)
    end

    return fig
end

# ╔═╡ 4cac6ee8-8c74-42f7-99f5-dab3804f74e3
vicsek_sum_weights()

# ╔═╡ Cell order:
# ╠═80b8650e-6390-11ef-2417-eb2b5b153f95
# ╠═c6939b67-42a4-45e1-9c61-33e9c1397d4c
# ╠═f788e3a3-181f-4676-9463-ed5b53e47223
# ╠═f68756e5-353a-4481-8e03-0796b5384f90
# ╠═00cf0514-7931-4d7e-a523-ae918eac6de3
# ╠═baa42609-1d76-439a-a82e-da277cbbc754
# ╠═ef2380b5-1ac8-43a7-be6c-b50d29608af5
# ╠═021d4f9e-5faa-4c22-9d2f-9ba0916d4b73
# ╠═7dbc6f7f-db53-45d6-8378-4dae562608f5
# ╠═f89e7cb5-f6f7-4fe0-adc0-96b10963e755
# ╠═a2e84401-9a34-4446-850b-74e729bc82b9
# ╠═4cac6ee8-8c74-42f7-99f5-dab3804f74e3
