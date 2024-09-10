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
    # const POINTTYPE = "equispaced"
    const POINTTYPE = "Chebyshev-1"
    # const POINTTYPE = "Chebyshev-2"
    # const POINTTYPE = "Gauss-Legendre"
    # const POINTTYPE = "Gauss-Lobatto"

    const MAXITER = 2048

    #! Data generation
    const SAVEDATA = false # take some time

    #! Ploting constants
    const ADDTITLE = true
    const SAVEPLOT = false
    const FONTSIZE = 20

    "Global parameters"
end

# ╔═╡ f788e3a3-181f-4676-9463-ed5b53e47223
function _save_weights(sas::src.SelfAffineSet, nb_pts_max::Int, suffix::String="")
    filename = sas.name * suffix

    open("../data-weights/$filename.csv", "w") do file
        write(file, "nb-points,sum-abs\n")
        write(file, "1,$(@sprintf("%.16e", 1))\n")

        for M in 2:nb_pts_max
            cbt = src.compute_cubature(sas, POINTTYPE, M; maxiter=MAXITER)

            n = @sprintf("%d", length(cbt))
            w = @sprintf("%.16e", sum(abs.(cbt.weights)))
            write(file, "$n,$w\n")
        end
    end

    return nothing
end

# ╔═╡ f68756e5-353a-4481-8e03-0796b5384f90
if SAVEDATA
    for (factor, suffix) in [(1 / 2, "1o2"), (1 / 3, "1o3"), (1 / 4, "1o4")]
        print((factor, suffix))
        _save_weights(src.cantor_set(factor, [0, 1]), 500, suffix)
        println(": done")
    end

    for (sas, suffix) in [
        (src.vicsek_2d(1 / 3), ""),
        (src.vicsek_2d(1 / 3, 0.4), "-0.4"),
        (src.vicsek_2d(1 / 3, π / 4), "-pio4"),
    ]
        _save_weights(sas, 32, suffix)
    end
end

# ╔═╡ 00cf0514-7931-4d7e-a523-ae918eac6de3
function cantor_set_weights(nb_pts::Int)
    sas = src.cantor_set(1 / 3, [0.0, 1.0])
    cbt = src.compute_cubature(sas, POINTTYPE, nb_pts; maxiter=MAXITER)
    w_lim = extrema(cbt.weights)
    W = maximum(abs.(w_lim))

    dx, dy = 0.025, 0.05 * (w_lim[2] - w_lim[1])
    X = (-dx, 1 + dx)
    Y = (w_lim[1] - dy, w_lim[2] + dy)

    pf = [src.vertices(sas.bounding_box)]
    for _ in 1:5
        pf = [S.(part) for S in sas.ifs for part in pf]
    end

    fig = Figure(; fontsize=FONTSIZE)

    ax_args::Dict{Symbol,Any} = Dict()
    if ADDTITLE
        ax_args[:title] = L"Cantor set with $\rho=1/3$ and $ M = %$nb_pts $"
        ax_args[:xlabel] = L"x"
        ax_args[:ylabel] = L"w"
    end
    ax = Axis(fig[1, 1]; ax_args...)

    limits!(ax, X, Y)

    for p in pf
        poly!(
            ax,
            [(p[1][1], Y[1]), (p[2][1], Y[1]), (p[2][1], Y[2]), (p[1][1], Y[2])];
            color=(:black, 0.3),
        )
    end

    scatter!(
        ax,
        reinterpret(Float64, cbt.points),
        cbt.weights;
        color=cbt.weights,
        colormap=:vik,
        colorrange=(-W, W),
        strokewidth=1,
        strokecolor=:black,
    )

    if SAVEPLOT
        save("1d-cantor-set-weights.pdf", fig)
    end

    return fig
end

# ╔═╡ baa42609-1d76-439a-a82e-da277cbbc754
cantor_set_weights(255)

# ╔═╡ ef2380b5-1ac8-43a7-be6c-b50d29608af5
function vicsek_weights(nb_pts::Int)
    sas = src.vicsek_2d(1 / 3, 0.4)
    cbt = src.compute_cubature(sas, POINTTYPE, nb_pts; maxiter=MAXITER)
    W = maximum(abs.(cbt.weights))

    pf = [[
        SVector{2,Float64}([-1, -1]),
        SVector{2,Float64}([1, -1]),
        SVector{2,Float64}([1, 1]),
        SVector{2,Float64}([-1, 1]),
    ]]
    for _ in 1:5
        pf = [S.(part) for S in sas.ifs for part in pf]
    end

    fig = Figure(; fontsize=FONTSIZE)

    ax_args::Dict{Symbol,Any} = Dict(:aspect => 1)
    if ADDTITLE
        M2 = nb_pts * nb_pts
        ax_args[:title] = L"Vicsek with $\rho=1/3$ and $ M = %$M2 $"
        ax_args[:xlabel] = L"x"
        ax_args[:ylabel] = L"y"
    end
    ax = Axis(fig[1, 1]; ax_args...)

    dx = 0.025
    limits!(ax, (-1 - dx, 1 + dx), (-1 - dx, 1 + dx))

    for part in pf
        poly!(ax, part; color=(:black, 0.5))
    end

    sc = scatter!(
        ax,
        cbt.points;
        color=cbt.weights,
        colormap=:vik,
        colorrange=(-W, W),
        strokewidth=1,
        strokecolor=:black,
    )

    if ADDTITLE
        Colorbar(fig[1, 2], sc; label=L"w")
    else
        Colorbar(fig[1, 2], sc)
    end

    if SAVEPLOT
        save("2d-vicsek-weights.pdf", fig)
    end

    return fig
end

# ╔═╡ 021d4f9e-5faa-4c22-9d2f-9ba0916d4b73
vicsek_weights(10)

# ╔═╡ 7dbc6f7f-db53-45d6-8378-4dae562608f5
function cantor_set_sum_weights()
    fig = Figure(; fontsize=FONTSIZE)

    ax_args::Dict{Symbol,Any} = Dict(:xlabel => L"M")
    if ADDTITLE
        ax_args[:title] = L"Sum of the absolute value of the weight for Cantor set$$"
        ax_args[:ylabel] = L"|\mathbf{w}|_1"
    end
    ax = Axis(fig[1, 1]; ax_args...)

    x = collect(1:500)
    y = 1 .+ 0.14 .* log.(x)
    lines!(ax, x, y; color=:black, label=L"\mathcal{O}(\log M)")

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

    if SAVEPLOT
        save("1d-cantor-set-sum-weights.pdf", fig)
    end

    return fig
end

# ╔═╡ f89e7cb5-f6f7-4fe0-adc0-96b10963e755
cantor_set_sum_weights()

# ╔═╡ a2e84401-9a34-4446-850b-74e729bc82b9
function vicsek_sum_weights()
    fig = Figure(; fontsize=FONTSIZE)

    ax_args::Dict{Symbol,Any} = Dict(:xlabel => L"M")
    if ADDTITLE
        ax_args[:title] = L"Sum of the absolute value of the weight for Vicsek$$"
        ax_args[:ylabel] = L"|\mathbf{w}|_1"
    end
    ax = Axis(fig[1, 1]; ax_args...)

    x = range(1, 1024; length=64)
    y = 1 .+ 0.007 .* log.(x) .^ 2
    lines!(ax, x, y; color=:black, label=L"\mathcal{O}(\log^2 M)")

    for (s, m, l) in [
        ("", :circle, L"\theta = 0"),
        ("-rot-0.4", :cross, L"\theta = 0.4"),
        ("-rot-pio4", :xcross, L"\theta = \pi / 4"),
    ]
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

        scatterlines!(ax, nb_pts, sum_abs; marker=m, linestyle=:dash, label=l)
    end

    axislegend(ax; position=:lt)

    if SAVEPLOT
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
