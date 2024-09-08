### A Pluto.jl notebook ###
# v0.19.46

using Markdown
using InteractiveUtils

# ╔═╡ aa732af4-6dc3-11ef-1ec5-1fed1fee0ea5
begin
    using Pkg: Pkg
    Pkg.activate(Base.current_project())
    Pkg.instantiate()

    using StaticArrays, CairoMakie

    import IFSCubature as src
end

# ╔═╡ 4b9f2377-8ba9-4c77-adbe-5d4f0dad22f6
begin
    #! Ploting constants
    const ADDTITLE = true
    const SAVEPLOT = false
    const FONTSIZE = 15

    "Global parameters"
end

# ╔═╡ d2aedcf8-321b-4747-812e-230d6117f9a6
function plot_nb_refine(
    sas::src.SelfAffineSet{3,T,9}, f0::src.Polyhedron{T}, nb_refine::Int, suffix::String=""
) where {T}
    fig = Figure()
    ax = Axis(fig[1, 1]; aspect=1)

    if SAVEPLOT
        save("$(sas.name)$suffix.pdf", fig)
    end

    return fig
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
    sas::src.SelfAffineSet{2,T,4}, f0::src.Polygon{T}, nb_refine::Int, suffix::String=""
) where {T}
    fig = Figure(; size=(600, 600), fontsize=FONTSIZE)

    ax_args::Dict{Symbol,Any} = Dict(:aspect => 1)
    if ADDTITLE
        ax_args[:title] = "$(sas.name)"
        ax_args[:xlabel] = L"x"
        ax_args[:ylabel] = L"y"
    end
    ax = Axis(fig[1, 1]; ax_args...)

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
    src.Polygon([[1.0, 0.0], [-0.5, √3 / 2], [-0.5, -√3 / 2]]),
    7,
)

# ╔═╡ fb4ab529-02ef-433e-9e39-db9d6a552e3b
plot_refine(
    src.koch_snowflake(),
    src.Polygon([[v for v in reverse(sincospi(2 * i//6))] for i in 0:5]),
    5,
)

# ╔═╡ 29cc3ac3-765a-47dd-96dc-27748d56cb53
plot_refine(
    src.vicsek_2d(1 / 3),
    src.Polygon([[1.0, 1.0], [-1.0, 1.0], [-1.0, -1.0], [1.0, -1.0]]),
    5,
)

# ╔═╡ 40af0eed-5d64-4b59-8a28-e1893812acc0
plot_refine(
    src.vicsek_2d(1 / 3, 0.4),
    src.Polygon([[1.0, 1.0], [-1.0, 1.0], [-1.0, -1.0], [1.0, -1.0]]),
    5,
    "-0.4",
)

# ╔═╡ 0ca96d6b-e308-4237-aa4d-325f4bdb9bcd
plot_refine(
    src.vicsek_2d(1 / 3, π / 4),
    src.Polygon([[1.0, 1.0], [-1.0, 1.0], [-1.0, -1.0], [1.0, -1.0]]),
    5,
    "-pio4",
)

# ╔═╡ b4f5039b-6f32-4d6f-ab93-7162d2989fbf
plot_refine(
    src.cantor_dust(1 / 3, [-1.0, 1.0], 2),
    src.Polygon([[1.0, 1.0], [-1.0, 1.0], [-1.0, -1.0], [1.0, -1.0]]),
    5,
)

# ╔═╡ d901e89b-c7f7-4dac-b815-bf4665bf5beb
begin
    local sas = src.cantor_dust_non_sym()
    local vs = src.fix_point.(sas.ifs)
    plot_refine(sas, src.Polygon(vs), 5)
end

# ╔═╡ 61e6c900-e44e-4b26-ab7c-d45a6c06e355
function plot_chaos_game(
    sas::src.SelfAffineSet{2,Float64,4}, nb_pts::Int, suffix::String=""
)
    fig = Figure(; size=(600, 600), fontsize=FONTSIZE)

    ax_args::Dict{Symbol,Any} = Dict(:aspect => 1)
    if ADDTITLE
        ax_args[:title] = "$(sas.name)"
        ax_args[:xlabel] = L"x"
        ax_args[:ylabel] = L"y"
    end
    ax = Axis(fig[1, 1]; ax_args...)

    box = sas.bounding_box
    r = max(box.paxis[1, 1], box.paxis[2, 2])
    _min = box.center .- r
    # _max = box.center .+ r

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
        (_min[1] .+ 0):((n - 1) .* h),
        (_min[2] .+ 0):((n - 1) .* h),
        f;
        colormap=Reverse(:grays),
        colorrange=(0, 1),
    )

    if SAVEPLOT
        save("$(sas.name)$suffix.png", fig)
    end

    return fig
end

# ╔═╡ a2919f15-267d-4e1a-bfa3-cfeac7751a54
plot_chaos_game(src.barnsley_fern(), 500_000)

# ╔═╡ Cell order:
# ╠═aa732af4-6dc3-11ef-1ec5-1fed1fee0ea5
# ╠═4b9f2377-8ba9-4c77-adbe-5d4f0dad22f6
# ╠═8ee12a1c-a9b4-4ab5-8e2f-e75e1aafa544
# ╠═e5da3ad5-0f2e-45f8-b79a-7fa765a924a2
# ╠═fb4ab529-02ef-433e-9e39-db9d6a552e3b
# ╠═a2919f15-267d-4e1a-bfa3-cfeac7751a54
# ╠═29cc3ac3-765a-47dd-96dc-27748d56cb53
# ╠═40af0eed-5d64-4b59-8a28-e1893812acc0
# ╠═0ca96d6b-e308-4237-aa4d-325f4bdb9bcd
# ╠═b4f5039b-6f32-4d6f-ab93-7162d2989fbf
# ╠═d901e89b-c7f7-4dac-b815-bf4665bf5beb
# ╠═d2aedcf8-321b-4747-812e-230d6117f9a6
# ╠═12794fa9-f6a7-457b-a44b-a902d0db00aa
# ╠═27c86a2c-e3af-4f66-8583-84302b042f40
# ╠═61e6c900-e44e-4b26-ab7c-d45a6c06e355
