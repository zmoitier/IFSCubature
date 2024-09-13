### A Pluto.jl notebook ###
# v0.19.46

using Markdown
using InteractiveUtils

# ╔═╡ 6d7e2666-71de-11ef-2330-1ff599d0d4b2
begin
    using Pkg: Pkg
    Pkg.activate(Base.current_project())
    Pkg.instantiate()

    using StaticArrays, CairoMakie

    import IFSCubature as src
end

# ╔═╡ 35bd0fbe-55f1-44d2-b9cf-b87903239e60
begin
    #! Type of points
    # const POINTTYPE = "equispaced"
    const POINTTYPE = "Chebyshev-1"
    # const POINTTYPE = "Chebyshev-2"
    # const POINTTYPE = "Gauss-Legendre"
    # const POINTTYPE = "Gauss-Lobatto"

    const FONTSIZE = 20

    "Global parameters"
end

# ╔═╡ 30c56128-b0e3-421e-a5d4-5c3eb83089c3
function plot_measure(sas::src.SelfAffineSet{1,T,1}, f_diam::Real, a::Real) where {T}
    hcbt = src.HCubature(src.barycenter_rule(sas), sas)
    while src.diameter(hcbt) > f_diam * src.diameter(sas)
        src.refine!(hcbt, sas)
    end

    fig = Figure(; size=(400, 400), fontsize=FONTSIZE)
    ax = Axis(fig[1, 1]; title=sas.name)

    ball = sas.bounding_ball
    for (_, S, μ) in hcbt.heap.valtree
        tmp = S(ball)
        lines!(ax, tmp.center[1] .+ [-tmp.radius, tmp.radius], [μ, μ]; color=(:black, a))
    end

    return fig
end

# ╔═╡ c926260d-a936-412d-82b9-9d8fe1e9b2f1
function plot_measure(sas::src.SelfAffineSet{2,T,4}, f_diam::Real) where {T}
    hcbt = src.HCubature(src.barycenter_rule(sas), sas)
    while src.diameter(hcbt) > f_diam * src.diameter(sas)
        src.refine!(hcbt, sas)
    end
    μ_max = maximum(μ for (_, _, μ) in hcbt.heap.valtree)

    fig = Figure(; size=(400, 400), fontsize=FONTSIZE)
    ax = Axis(fig[1, 1]; aspect=1, title=sas.name)

    ball = sas.bounding_ball
    while !isempty(hcbt.heap)
        _, S, μ = pop!(hcbt.heap)
        tmp = S(ball)
        poly!(
            ax,
            Circle(Point2f(tmp.center), tmp.radius * 0.75);
            color=μ,
            colormap=:matter,
            colorrange=(0, μ_max),
        )
    end
    Colorbar(fig[1, 2]; limits=(0, μ_max), colormap=:matter)

    return fig
end

# ╔═╡ 46ad88ba-a14d-4164-bfcd-8910c3fd5ddd
plot_measure(src.cantor_set([0.5, 0.5], [0.0, 1.0], [0.75, 0.25]), 0.01, 1)

# ╔═╡ ba51ee12-f6ec-48e3-be4c-89f5923c1177
plot_measure(src.cantor_dust(1 / 3, [0.0, 1.0], 2), 0.05)

# ╔═╡ 9b37f854-b180-4346-8a46-6a3c47ee8379
plot_measure(src.sierpinski_triangle(), 0.05)

# ╔═╡ 18194164-9b4b-4b31-8221-fdf94fe2ac02
plot_measure(src.sierpinski_triangle_fat(2), 0.05)

# ╔═╡ c675f766-0400-40d5-ad88-37c8068be772
plot_measure(src.vicsek_2d(1 / 3), 0.05)

# ╔═╡ 44d931e8-e1f2-49b5-9b8d-e0a317fe19ab
plot_measure(src.sierpinski_carpet(), 0.05)

# ╔═╡ 03d9bc59-b5eb-4463-a2fc-6770cb14cfc5
plot_measure(src.koch_snowflake(), 0.1)

# ╔═╡ 5e7fc012-c509-4235-8521-ce512ba15116
plot_measure(src.gosper_flowsnake(), 0.1)

# ╔═╡ 2a031615-9dc7-446c-a141-ea2fa619a370
plot_measure(src.brick_2d(), 0.1)

# ╔═╡ 9faa03aa-6158-48f8-828c-bfceac5ef341
plot_measure(src.durer_pentagon(), 0.1)

# ╔═╡ ac36e0d5-2a0b-4740-a50e-d808bc9c1f58
plot_measure(src.fudgeflake(), 0.1)

# ╔═╡ 75dc6b84-8044-4622-8708-df40069aebd2
plot_measure(src.heighway_dragon(), 0.1)

# ╔═╡ 3c35046c-9135-4e56-a002-eed33b96773c
plot_measure(src.levy_dragon(), 0.1)

# ╔═╡ c92f63b9-9e71-4ee9-90b3-8aa643d29ba9
plot_measure(src.terdragon(), 0.1)

# ╔═╡ 1349443e-2e12-4941-8dc5-0a0d54b99682
plot_measure(src.twindragon(), 0.1)

# ╔═╡ ff1672f3-870f-4e5a-aa98-bb4a03dca706
plot_measure(src.cantor_dust_non_sym(), 0.05)

# ╔═╡ a6e4a459-084a-4b44-b450-d5b5bd4c1e06
plot_measure(src.barnsley_fern(), 0.05)

# ╔═╡ Cell order:
# ╠═6d7e2666-71de-11ef-2330-1ff599d0d4b2
# ╠═35bd0fbe-55f1-44d2-b9cf-b87903239e60
# ╠═46ad88ba-a14d-4164-bfcd-8910c3fd5ddd
# ╠═ba51ee12-f6ec-48e3-be4c-89f5923c1177
# ╠═9b37f854-b180-4346-8a46-6a3c47ee8379
# ╠═18194164-9b4b-4b31-8221-fdf94fe2ac02
# ╠═c675f766-0400-40d5-ad88-37c8068be772
# ╠═44d931e8-e1f2-49b5-9b8d-e0a317fe19ab
# ╠═03d9bc59-b5eb-4463-a2fc-6770cb14cfc5
# ╠═5e7fc012-c509-4235-8521-ce512ba15116
# ╠═2a031615-9dc7-446c-a141-ea2fa619a370
# ╠═9faa03aa-6158-48f8-828c-bfceac5ef341
# ╠═ac36e0d5-2a0b-4740-a50e-d808bc9c1f58
# ╠═75dc6b84-8044-4622-8708-df40069aebd2
# ╠═3c35046c-9135-4e56-a002-eed33b96773c
# ╠═c92f63b9-9e71-4ee9-90b3-8aa643d29ba9
# ╠═1349443e-2e12-4941-8dc5-0a0d54b99682
# ╠═ff1672f3-870f-4e5a-aa98-bb4a03dca706
# ╠═a6e4a459-084a-4b44-b450-d5b5bd4c1e06
# ╠═30c56128-b0e3-421e-a5d4-5c3eb83089c3
# ╠═c926260d-a936-412d-82b9-9d8fe1e9b2f1
