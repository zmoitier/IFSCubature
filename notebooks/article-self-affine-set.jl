### A Pluto.jl notebook ###
# v0.20.21

using Markdown
using InteractiveUtils

# ╔═╡ aa732af4-6dc3-11ef-1ec5-1fed1fee0ea5
begin
    using Pkg: Pkg
    Pkg.activate(Base.current_project())
    Pkg.instantiate()

    using LinearAlgebra, StaticArrays, CairoMakie

    import IFSCubature as src
end

# ╔═╡ 4b9f2377-8ba9-4c77-adbe-5d4f0dad22f6
begin
    #! Ploting constants
    const ADDTITLE = false
    const SAVEPLOT = true
    const FONTSIZE = 15

    "Global parameters"
end

# ╔═╡ 12794fa9-f6a7-457b-a44b-a902d0db00aa
function plot_refine(
    sas::src.SelfAffineSet{1,T,1}, f0::src.Segment{T}, nb_refine::Int, suffix::String=""
) where {T}
    fig = Figure(; fontsize=FONTSIZE)

    ax_args::Dict{Symbol,Any} = Dict(:aspect => 1)
    if ADDTITLE
        ax_args[:title] = "$(sas.name)"
        ax_args[:xlabel] = L"x"
        ax_args[:ylabel] = L"p"
    end
    ax = Axis(fig[1, 1]; ax_args...)

    fp = [[f0.vertices[1], f0.vertices[2]]]
    args = Dict(:color => 1, :colormap => :tab10, :colorrange => (1, 10), :linewidth => 3)

    lines!(ax, fp[1], fill(0, 2); args...)
    for p in 1:nb_refine
        fp = [S.(part) for S in sas.ifs for part in fp]
        for part in fp
            lines!(ax, part, fill(-p, 2); args...)
        end
    end

    if SAVEPLOT
        save("$(sas.name)$suffix.pdf", fig)
    end

    return fig
end

# ╔═╡ 27c86a2c-e3af-4f66-8583-84302b042f40
function plot_refine(
    sas::src.SelfAffineSet{2,T,4},
    f0::src.Polygon{T};
    nb_refine::Int=1,
    α::Real=0.5,
    size::Tuple{Int,Int}=(600, 600),
    suffix::String="",
) where {T}
    box = sas.bounding_box
    r = diag(box.paxis)
    _min = box.center .- r
    _max = box.center .+ r

    a, b = _max - _min
    δ = 0.025 * norm(_max - _min)

    fig = Figure(; size=size, fontsize=FONTSIZE, figure_padding=1)

    ax_args::Dict{Symbol,Any} = Dict(:aspect => a / b)
    if ADDTITLE
        ax_args[:title] = "$(sas.name)"
        ax_args[:xlabel] = L"x"
        ax_args[:ylabel] = L"y"
    end
    ax = Axis(fig[1, 1]; ax_args...)

    xlims!(ax, (_min[1] - δ, _max[1] + δ))
    ylims!(ax, (_min[2] - δ, _max[2] + δ))

    poly!(ax, Point2f.(src.vertices(box))[[1, 2, 4, 3]]; color=(:black, α))

    fp = [f0.vertices]
    for _ in 1:nb_refine
        fp = [S.(pts) for S in sas.ifs for pts in fp]
    end

    for pts in fp
        poly!(ax, pts; color=:black)
    end

    if SAVEPLOT
        save("$(sas.name)$suffix.pdf", fig)
    end

    return fig
end

# ╔═╡ 8ee12a1c-a9b4-4ab5-8e2f-e75e1aafa544
plot_refine(src.cantor_set(1 / 3, [0.0, 1.0]), src.Segment(0.0, 1.0), 4)

# ╔═╡ e5da3ad5-0f2e-45f8-b79a-7fa765a924a2
plot_refine(
    src.sierpinski_triangle_fat(2),
    src.Polygon([[1.0, 0.0], [-0.5, √3 / 2], [-0.5, -√3 / 2]]);
    nb_refine=7,
    α=0.1,
    size=(550, 600),
)

# ╔═╡ 29cc3ac3-765a-47dd-96dc-27748d56cb53
plot_refine(
    src.vicsek_2d(1 / 3),
    src.Polygon([[1.0, 1.0], [-1.0, 1.0], [-1.0, -1.0], [1.0, -1.0]]);
    nb_refine=5,
    α=0.05,
    size=(625, 600),
)

# ╔═╡ 40af0eed-5d64-4b59-8a28-e1893812acc0
plot_refine(
    src.vicsek_2d(1 / 3, 0.4),
    src.Polygon([[1.0, 1.0], [-1.0, 1.0], [-1.0, -1.0], [1.0, -1.0]]);
    nb_refine=5,
    α=0.05,
    size=(625, 600),
    suffix="-0.4",
)

# ╔═╡ 0ca96d6b-e308-4237-aa4d-325f4bdb9bcd
plot_refine(
    src.vicsek_2d(1 / 3, π / 4),
    src.Polygon([[1.0, 1.0], [-1.0, 1.0], [-1.0, -1.0], [1.0, -1.0]]);
    nb_refine=5,
    α=0.05,
    size=(625, 600),
    suffix="-pio4",
)

# ╔═╡ b4f5039b-6f32-4d6f-ab93-7162d2989fbf
plot_refine(
    src.cantor_dust(1 / 3, [-1.0, 1.0], 2),
    src.Polygon([[1.0, 1.0], [-1.0, 1.0], [-1.0, -1.0], [1.0, -1.0]]);
    nb_refine=5,
    α=0.05,
    size=(625, 600),
)

# ╔═╡ d901e89b-c7f7-4dac-b815-bf4665bf5beb
begin
    local sas = src.cantor_dust_non_sym()
    local vs = src.fix_point.(sas.ifs)
    plot_refine(sas, src.Polygon(vs); nb_refine=5, α=0.05, size=(750, 600))
end

# ╔═╡ 61e6c900-e44e-4b26-ab7c-d45a6c06e355
function plot_chaos_game(
    sas::src.SelfAffineSet{2,Float64,4};
    nb_pts::Int=1_024,
    α::Real=0.5,
    size::Tuple{Int,Int}=(600, 600),
    suffix::String="",
)
    box = sas.bounding_box
    r = diag(box.paxis)
    _min = box.center .- r
    _max = box.center .+ r

    a, b = _max - _min
    δ = 0.025 * norm(_max - _min)

    fig = Figure(; size=size, fontsize=FONTSIZE, figure_padding=1)

    ax_args::Dict{Symbol,Any} = Dict(:aspect => a / b)
    if ADDTITLE
        ax_args[:title] = "$(sas.name)"
        ax_args[:xlabel] = L"x"
        ax_args[:ylabel] = L"y"
    end
    ax = Axis(fig[1, 1]; ax_args...)

    xlims!(ax, (_min[1] - δ, _max[1] + δ))
    ylims!(ax, (_min[2] - δ, _max[2] + δ))

    poly!(ax, Point2f.(src.vertices(box))[[1, 2, 4, 3]]; color=(:black, α))

    n = 512
    h = 2 * r / (n - 1)

    function xy_to_ij(xy)
        return floor.(Int, (xy - _min) ./ h .+ 0.5) .+ 1
    end

    p = cumsum(sas.measure)

    f = fill(NaN, n, n)
    for c in src.fix_points(sas)
        xy = MVector{2,Float64}(c)
        ij = xy_to_ij(xy)
        f[ij[1], ij[2]] = 1
        for _ in 1:(nb_pts ÷ length(sas.ifs))
            k = searchsortedfirst(p, rand())
            xy = sas.ifs[k](xy)
            ij = xy_to_ij(xy)
            f[ij[1], ij[2]] = 1
        end
    end

    heatmap!(
        ax,
        range(_min[1], _max[1], n),
        range(_min[2], _max[2], n),
        f;
        colormap=Reverse(:grays),
        colorrange=(0, 1),
    )

    if SAVEPLOT
        save("$(sas.name)$suffix.png", fig)
    end

    return fig
end

# ╔═╡ e61bb2e3-0a03-4732-ba8d-3845b5439768
# plot_refine(
#     src.koch_snowflake(),
#     src.Polygon([[v for v in reverse(sincospi(2 * i//6))] for i in 0:5]),
#     nb_refine=5,
#     α=0.1,
# 	size=(725,600),
# )
plot_chaos_game(src.koch_snowflake(); nb_pts=2_000_000, α=0.1, size=(725, 600))

# ╔═╡ 3c4d8a45-2cc7-4f5c-b842-6ee6ae09da95
plot_chaos_game(src.barnsley_fern(); nb_pts=2_000_000, α=0.1, size=(350, 600))

# ╔═╡ Cell order:
# ╠═aa732af4-6dc3-11ef-1ec5-1fed1fee0ea5
# ╠═4b9f2377-8ba9-4c77-adbe-5d4f0dad22f6
# ╠═8ee12a1c-a9b4-4ab5-8e2f-e75e1aafa544
# ╠═e5da3ad5-0f2e-45f8-b79a-7fa765a924a2
# ╠═29cc3ac3-765a-47dd-96dc-27748d56cb53
# ╠═40af0eed-5d64-4b59-8a28-e1893812acc0
# ╠═0ca96d6b-e308-4237-aa4d-325f4bdb9bcd
# ╠═b4f5039b-6f32-4d6f-ab93-7162d2989fbf
# ╠═d901e89b-c7f7-4dac-b815-bf4665bf5beb
# ╠═e61bb2e3-0a03-4732-ba8d-3845b5439768
# ╠═3c4d8a45-2cc7-4f5c-b842-6ee6ae09da95
# ╠═12794fa9-f6a7-457b-a44b-a902d0db00aa
# ╠═27c86a2c-e3af-4f66-8583-84302b042f40
# ╠═61e6c900-e44e-4b26-ab7c-d45a6c06e355
