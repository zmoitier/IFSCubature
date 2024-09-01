### A Pluto.jl notebook ###
# v0.19.46

using Markdown
using InteractiveUtils

# ╔═╡ 09475398-6441-11ef-3d59-05d1f8d5d356
begin
    using Pkg: Pkg
    Pkg.activate(Base.current_project())
    Pkg.instantiate()

    using LinearAlgebra, CairoMakie

    import IFSCubature as src
end

# ╔═╡ a19083f3-fc05-42eb-bda3-7ae44a98e382
function _plot(attractor::src.SelfAffineSet{3,T,9}, p_max::Int) where {T} end

# ╔═╡ c241a629-26bf-450d-a121-87e7d649fd94
function _plot(sas::src.SelfAffineSet{1,T,1}, f0::src.Segment{T}, p_max::Int) where {T}
    fig = Figure()
    ax = Axis(fig[1, 1]; xlabel=L"x", aspect=1)

    box = [x[1] for x in src.vertices(sas.bounding_box)]
    for p in 0:p_max
        lines!(ax, box, fill(-p, 2); linewidth=16, color=(:black, 0.25))
    end

    fp = [[f0.vertices[1], f0.vertices[2]]]
    args = Dict(:color => 1, :colormap => :tab10, :colorrange => (1, 10))
    scatterlines!(ax, fp[1], fill(0, 2); args...)
    for p in 1:p_max
        fp = [S.(part) for S in sas.ifs for part in fp]
        for part in fp
            scatterlines!(ax, part, fill(-p, 2); args...)
        end
    end

    return fig
end

# ╔═╡ 747e7ba7-ff26-48f0-861b-3a286dfc1d65
function _plot(sas::src.SelfAffineSet{2,T,4}, f0::src.Polygon{T}, p_max::Int) where {T}
    tab10 = Makie.to_colormap(:tab10)

    fig = Figure(; size=(400, p_max * 500))
    axs = [
        Axis(fig[p, 1]; xlabel=L"x", ylabel=L"y", aspect=1, title=L"p = %$p") for
        p in 0:p_max
    ]

    ball = sas.bounding_ball
    box = Point2f.(src.vertices(sas.bounding_box))[[1, 2, 4, 3]]

    fp = [f0.vertices]
    poly!(axs[1], f0.vertices; color=(:black, 0.25), strokecolor=:black, strokewidth=2)
    for ax in axs[2:end]
        fp = [S.(part) for S in sas.ifs for part in fp]
        for part in fp
            poly!(ax, part; color=(:black, 0.25), strokecolor=:black, strokewidth=2)
        end
    end

    for ax in axs
        poly!(
            ax,
            Circle(Point2f(ball.center), ball.radius);
            color=(:black, 0),
            strokecolor=(tab10[1], 0.75),
            strokewidth=2,
        )
        poly!(ax, box; color=(:black, 0), strokecolor=(tab10[2], 0.75), strokewidth=2)
    end

    return fig
end

# ╔═╡ 05935561-cc1c-43fa-ae37-4e9c447217a0
_plot(src.cantor_set(1 / 3, [0.0, 1.0]), src.Segment(0.0, 1.0), 3)

# ╔═╡ ea427ad6-71a8-45f6-bd58-4aa3c8c90fb0
_plot(
    src.cantor_dust(1 / 3, [-1.0, 1.0], 2),
    src.Polygon([[1.0, 1.0], [-1.0, 1.0], [-1.0, -1.0], [1.0, -1.0]]),
    2,
)

# ╔═╡ d22e02bd-6bec-48cb-a3e0-236ef67bd0dd
_plot(
    src.sierpinski_triangle(), src.Polygon([[1.0, 0.0], [-0.5, √3 / 2], [-0.5, -√3 / 2]]), 2
)

# ╔═╡ aeca916d-a93d-4aaf-9dd2-e2e7e6c54c93
_plot(
    src.fat_sierpinski_triangle(2),
    src.Polygon([[1.0, 0.0], [-0.5, √3 / 2], [-0.5, -√3 / 2]]),
    2,
)

# ╔═╡ 23c0f5f1-2319-4da8-8877-16515e01ad7f
_plot(
    src.vicsek_2d(1 / 3),
    src.Polygon([[1.0, 1.0], [-1.0, 1.0], [-1.0, -1.0], [1.0, -1.0]]),
    2,
)

# ╔═╡ 383424fd-98f8-42df-9337-b4d9feb9ce64
_plot(
    src.vicsek_2d(1 / 3, π / 4),
    src.Polygon([[1.0, 1.0], [-1.0, 1.0], [-1.0, -1.0], [1.0, -1.0]]),
    2,
)

# ╔═╡ 9ad71f58-6ed0-4fe0-9974-30a2ba4c3e8f
_plot(
    src.sierpinski_carpet(),
    src.Polygon([[1.0, 1.0], [-1.0, 1.0], [-1.0, -1.0], [1.0, -1.0]]),
    2,
)

# ╔═╡ 7286d0a4-5d33-40ce-bd00-8a39f10aac25
_plot(
    src.koch_snowflake(),
    src.Polygon([[v for v in reverse(sincospi(2 * i//6))] for i in 0:5]),
    2,
)

# ╔═╡ def9d8ca-ebad-4d32-bdfe-7c09e9ceaef4
_plot(
    src.gosper_flowsnake(),
    src.Polygon([[v for v in reverse(sincospi(2 * i//6))] for i in 0:5]),
    2,
)

# ╔═╡ 4ec105bf-e55e-4581-8f12-a924f9855360
begin
    local sas = src.cantor_dust_non_sym()
    local vs = src.fix_point.(sas.ifs)
    _plot(sas, src.Polygon(vs), 2)
end

# ╔═╡ 2347a278-0439-416f-b5a2-8d1ff151e224
begin
    local sas = src.barnsley_fern()
    local vs = src.fix_point.(sas.ifs)
    _plot(sas, src.Polygon(vs), 3)
end

# ╔═╡ Cell order:
# ╠═09475398-6441-11ef-3d59-05d1f8d5d356
# ╠═05935561-cc1c-43fa-ae37-4e9c447217a0
# ╠═ea427ad6-71a8-45f6-bd58-4aa3c8c90fb0
# ╠═d22e02bd-6bec-48cb-a3e0-236ef67bd0dd
# ╠═aeca916d-a93d-4aaf-9dd2-e2e7e6c54c93
# ╠═23c0f5f1-2319-4da8-8877-16515e01ad7f
# ╠═383424fd-98f8-42df-9337-b4d9feb9ce64
# ╠═9ad71f58-6ed0-4fe0-9974-30a2ba4c3e8f
# ╠═7286d0a4-5d33-40ce-bd00-8a39f10aac25
# ╠═def9d8ca-ebad-4d32-bdfe-7c09e9ceaef4
# ╠═4ec105bf-e55e-4581-8f12-a924f9855360
# ╠═2347a278-0439-416f-b5a2-8d1ff151e224
# ╠═a19083f3-fc05-42eb-bda3-7ae44a98e382
# ╠═c241a629-26bf-450d-a121-87e7d649fd94
# ╠═747e7ba7-ff26-48f0-861b-3a286dfc1d65
