### A Pluto.jl notebook ###
# v0.20.20

using Markdown
using InteractiveUtils

# ╔═╡ 49a64b38-6888-11ef-05ba-23348a353dd7
begin
    using Pkg: Pkg
    Pkg.activate(Base.current_project())
    Pkg.instantiate()

    using StaticArrays, CairoMakie

    import IFSCubature as src
end

# ╔═╡ ae4553bc-6cc2-4595-805f-c42dcd1b8573
function plot_bounding_3d(;
    sas::src.SelfAffineSet{3,T,9}, p_max::Int=0, nb_pts_chaos=10_000, α_attractor::Real=0.75
) where {T} end

# ╔═╡ 5fd68c56-a139-4fb3-952a-a42b05f15f33
plot_bounding_3d(; sas=src.cantor_dust(1 / 3, [-1.0, 1.0], 3), p_max=0, α_attractor=0.75)

# ╔═╡ aef7e47f-f688-4dc3-ad37-b028a7084dc8
plot_bounding_3d(; sas=src.sierpinski_tetrahedron(), p_max=0, α_attractor=0.75)

# ╔═╡ 099f41c1-0a10-4c8c-af5e-d420dd71a890
plot_bounding_3d(; sas=src.menger_sponge(), p_max=0, α_attractor=0.75)

# ╔═╡ ea213a70-0d1e-47f7-b5d3-b80f197d4696
plot_bounding_3d(; sas=src.vicsek_3d(1 / 3, true), p_max=0, α_attractor=0.75)

# ╔═╡ 5106bc48-4ae3-43da-857b-d64fc553a351
plot_bounding_3d(; sas=src.brick_3d(), p_max=0, α_attractor=0.75)

# ╔═╡ 787f1d23-9b60-4111-81c0-e212d492ba10
function plot_chaos_game!(
    ax::Axis, sas::src.SelfAffineSet{D,T,N}, nb_pts::Int
) where {D,T,N}
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
        xy = MVector{D,T}(c)
        ij = xy_to_ij(xy)
        f[ij[1], ij[2]] = 1
        for _ in 1:(nb_pts ÷ length(sas.ifs))
            k = searchsortedfirst(p, rand())
            xy = sas.ifs[k](xy)
            ij = xy_to_ij(xy)
            f[ij[1], ij[2]] = 1
        end
    end

    v = collect(0:(n - 1)) .* h
    heatmap!(ax, _min[1] .+ v, _min[2] .+ v, f; colormap=Reverse(:grays), colorrange=(0, 1))

    return nothing
end

# ╔═╡ e77088ae-f27b-4b6c-bd1c-881068a95e3b
function plot_bounding_2d(;
    sas::src.SelfAffineSet{2,T,4}, p_max::Int=0, nb_pts_chaos=10_000, α_attractor::Real=0.75
) where {T}
    colors = Makie.wong_colors()

    fig = Figure()
    ax = Axis(fig[1, 1]; xlabel=L"x", ylabel=L"y", aspect=DataAspect(), title=sas.name)

    plot_chaos_game!(ax, sas, nb_pts_chaos)

    balls = [sas.bounding_ball]
    ball_style = Dict(
        :color => (:black, 0), :strokecolor => (colors[1], 0.5), :strokewidth => 2
    )

    boxes = [sas.bounding_box]
    box_style = Dict(
        :color => (:black, 0), :strokecolor => (colors[2], 0.5), :strokewidth => 2
    )

    poly!(ax, Circle(Point2f(balls[1].center), balls[1].radius); ball_style...)
    poly!(ax, Point2f.(src.vertices(boxes[1]))[[1, 2, 4, 3]]; box_style...)

    for _ in 1:p_max
        balls = [S(ball) for S in sas.ifs for ball in balls]
        for ball in balls
            poly!(ax, Circle(Point2f(ball.center), ball.radius); ball_style...)
        end
        boxes = [S(box) for S in sas.ifs for box in boxes]
        for box in boxes
            poly!(ax, Point2f.(src.vertices(box))[[1, 2, 4, 3]]; box_style...)
        end
    end

    return fig
end

# ╔═╡ 3f57b6f2-9e65-4ede-a0ef-c4b042831e60
plot_bounding_2d(; sas=src.cantor_dust(1 / 3, [-1.0, 1.0], 2), p_max=0)

# ╔═╡ 35b630b2-ee17-4708-a0b5-9a68e2137fa8
plot_bounding_2d(; sas=src.sierpinski_triangle(), p_max=0, α_attractor=0.75)

# ╔═╡ 6a8f45fd-d2b7-492e-a16c-62b0c8422a1d
plot_bounding_2d(; sas=src.sierpinski_triangle_fat(2), p_max=0, α_attractor=0.75)

# ╔═╡ 2280dbd2-889a-4a38-8c98-e622dc7ab53d
plot_bounding_2d(; sas=src.vicsek_2d(1 / 3), p_max=0, α_attractor=0.75)

# ╔═╡ f6788e2d-ad6f-4dfc-9676-420c3a94f84e
plot_bounding_2d(; sas=src.vicsek_2d(1 / 3, 0.4), p_max=0, α_attractor=0.75)

# ╔═╡ ac5bf72a-9274-4ea0-bf25-3444a7367364
plot_bounding_2d(; sas=src.sierpinski_carpet(), p_max=0, α_attractor=0.75)

# ╔═╡ 23bb9222-c6e5-41fd-8aa6-1ba55e748c3e
plot_bounding_2d(; sas=src.koch_snowflake(), p_max=0, α_attractor=0.75)

# ╔═╡ 2a671d8a-c2c6-4bf1-a84b-ff690ba52454
plot_bounding_2d(; sas=src.gosper_flowsnake(), p_max=0, α_attractor=0.75)

# ╔═╡ 2757e2a6-16d2-4083-9166-a553da04c63d
plot_bounding_2d(; sas=src.brick_2d(), p_max=0, α_attractor=0.75)

# ╔═╡ 35c82cbd-4622-431c-99fc-8b90e782f46f
plot_bounding_2d(; sas=src.durer_pentagon(), p_max=0, α_attractor=0.75)

# ╔═╡ bdf97855-7344-4bb0-ac5e-5bd0a98aa3db
plot_bounding_2d(; sas=src.fudgeflake(), p_max=0, α_attractor=0.75)

# ╔═╡ 7e75dd41-c44b-4fef-9d11-f20549a4bb3e
plot_bounding_2d(; sas=src.heighway_dragon(), p_max=0, α_attractor=0.75)

# ╔═╡ c657bba6-4366-496a-bdf2-12209f0c6214
plot_bounding_2d(; sas=src.levy_dragon(), p_max=0, α_attractor=0.75)

# ╔═╡ e3b08dee-9180-4dcc-9014-739de1358fbc
plot_bounding_2d(; sas=src.terdragon(), p_max=0, α_attractor=0.75)

# ╔═╡ 428a956e-cf29-4b8b-a812-6f6cd89a91c7
plot_bounding_2d(; sas=src.twindragon(), p_max=0, α_attractor=0.75)

# ╔═╡ fdc89186-6d96-488b-8bf2-cf97ae18a54a
plot_bounding_2d(; sas=src.cantor_dust_non_sym(), p_max=0, α_attractor=0.75)

# ╔═╡ 80907e2f-a22a-4a6d-9fb6-221cecc6fca5
plot_bounding_2d(; sas=src.barnsley_fern(), p_max=0, α_attractor=0.75)

# ╔═╡ Cell order:
# ╠═49a64b38-6888-11ef-05ba-23348a353dd7
# ╠═3f57b6f2-9e65-4ede-a0ef-c4b042831e60
# ╠═35b630b2-ee17-4708-a0b5-9a68e2137fa8
# ╠═6a8f45fd-d2b7-492e-a16c-62b0c8422a1d
# ╠═2280dbd2-889a-4a38-8c98-e622dc7ab53d
# ╠═f6788e2d-ad6f-4dfc-9676-420c3a94f84e
# ╠═ac5bf72a-9274-4ea0-bf25-3444a7367364
# ╠═23bb9222-c6e5-41fd-8aa6-1ba55e748c3e
# ╠═2a671d8a-c2c6-4bf1-a84b-ff690ba52454
# ╠═2757e2a6-16d2-4083-9166-a553da04c63d
# ╠═35c82cbd-4622-431c-99fc-8b90e782f46f
# ╠═bdf97855-7344-4bb0-ac5e-5bd0a98aa3db
# ╠═7e75dd41-c44b-4fef-9d11-f20549a4bb3e
# ╠═c657bba6-4366-496a-bdf2-12209f0c6214
# ╠═e3b08dee-9180-4dcc-9014-739de1358fbc
# ╠═428a956e-cf29-4b8b-a812-6f6cd89a91c7
# ╠═fdc89186-6d96-488b-8bf2-cf97ae18a54a
# ╠═80907e2f-a22a-4a6d-9fb6-221cecc6fca5
# ╠═5fd68c56-a139-4fb3-952a-a42b05f15f33
# ╠═aef7e47f-f688-4dc3-ad37-b028a7084dc8
# ╠═099f41c1-0a10-4c8c-af5e-d420dd71a890
# ╠═ea213a70-0d1e-47f7-b5d3-b80f197d4696
# ╠═5106bc48-4ae3-43da-857b-d64fc553a351
# ╠═e77088ae-f27b-4b6c-bd1c-881068a95e3b
# ╠═ae4553bc-6cc2-4595-805f-c42dcd1b8573
# ╠═787f1d23-9b60-4111-81c0-e212d492ba10
