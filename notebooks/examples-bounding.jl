### A Pluto.jl notebook ###
# v0.19.46

using Markdown
using InteractiveUtils

# ╔═╡ 49a64b38-6888-11ef-05ba-23348a353dd7
begin
    using Pkg: Pkg
    Pkg.activate(Base.current_project())
    Pkg.instantiate()

    using CairoMakie

    import IFSCubature as src
end

# ╔═╡ 1005ccc7-156f-4ea1-a60a-46497e8f144c
begin
    const p_max = 1

    "Global parameters"
end

# ╔═╡ 912bc5b6-a174-4ed0-9fa1-5f4e8de4d5fb
function _plot(sas::src.SelfAffineSet{3,T,9}, p_max::Int) where {T} end

# ╔═╡ 7436e59c-1b01-46cf-b6c6-0add22c4b7d4
function _plot(sas::src.SelfAffineSet{1,T,1}, p_max::Int) where {T}
    tab10 = Makie.to_colormap(:tab10)

    fig = Figure()
    ax = Axis(fig[1, 1]; xlabel=L"x", aspect=1, title=sas.name)

    balls = [sas.bounding_ball]
    ball_style = Dict(:color => tab10[1])

    boxes = [sas.bounding_box]
    box_style = Dict(:color => tab10[2], :linestyle => :dash)

    lines!(ax, src.vertices(balls[1]), fill(0, 2); ball_style...)
    lines!(ax, src.vertices(boxes[1]), fill(0, 2); box_style...)

    for p in 1:p_max
        balls = [S(ball) for S in sas.ifs for ball in balls]
        for ball in balls
            lines!(ax, src.vertices(ball), fill(-p, 2); ball_style...)
        end
        boxes = [S(box) for S in sas.ifs for box in boxes]
        for box in boxes
            lines!(ax, src.vertices(box), fill(-p, 2); box_style...)
        end
    end

    return fig
end

# ╔═╡ e77088ae-f27b-4b6c-bd1c-881068a95e3b
function _plot(sas::src.SelfAffineSet{2,T,4}, p_max::Int) where {T}
    tab10 = Makie.to_colormap(:tab10)

    fig = Figure()
    ax = Axis(fig[1, 1]; xlabel=L"x", ylabel=L"y", aspect=1, title=sas.name)

    balls = [sas.bounding_ball]
    ball_style = Dict(
        :color => (:black, 0), :strokecolor => (tab10[1], 0.5), :strokewidth => 2
    )

    boxes = [sas.bounding_box]
    box_style = Dict(
        :color => (:black, 0), :strokecolor => (tab10[2], 0.5), :strokewidth => 2
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

# ╔═╡ 50acd7c4-0770-4b11-aea4-f852236e4f31
_plot(src.cantor_set(1 / 3, [0.0, 1.0]), p_max)

# ╔═╡ 3f57b6f2-9e65-4ede-a0ef-c4b042831e60
_plot(src.cantor_dust(1 / 3, [-1.0, 1.0], 2), p_max)

# ╔═╡ 35b630b2-ee17-4708-a0b5-9a68e2137fa8
_plot(src.sierpinski_triangle(), p_max)

# ╔═╡ 6a8f45fd-d2b7-492e-a16c-62b0c8422a1d
_plot(src.sierpinski_triangle_fat(2), p_max)

# ╔═╡ 2280dbd2-889a-4a38-8c98-e622dc7ab53d
_plot(src.vicsek_2d(1 / 3), p_max)

# ╔═╡ f6788e2d-ad6f-4dfc-9676-420c3a94f84e
_plot(src.vicsek_2d(1 / 3, π / 4), p_max)

# ╔═╡ ac5bf72a-9274-4ea0-bf25-3444a7367364
_plot(src.sierpinski_carpet(), p_max)

# ╔═╡ 23bb9222-c6e5-41fd-8aa6-1ba55e748c3e
_plot(src.koch_snowflake(), p_max)

# ╔═╡ 2a671d8a-c2c6-4bf1-a84b-ff690ba52454
_plot(src.gosper_flowsnake(), p_max)

# ╔═╡ 2757e2a6-16d2-4083-9166-a553da04c63d
_plot(src.brick_2d(), p_max)

# ╔═╡ 35c82cbd-4622-431c-99fc-8b90e782f46f
_plot(src.durer_pentagon(), p_max)

# ╔═╡ bdf97855-7344-4bb0-ac5e-5bd0a98aa3db
_plot(src.fudgeflake(), p_max)

# ╔═╡ 7e75dd41-c44b-4fef-9d11-f20549a4bb3e
_plot(src.heighway_dragon(), p_max)

# ╔═╡ c657bba6-4366-496a-bdf2-12209f0c6214
_plot(src.levy_dragon(), p_max)

# ╔═╡ e3b08dee-9180-4dcc-9014-739de1358fbc
_plot(src.terdragon(), p_max)

# ╔═╡ 428a956e-cf29-4b8b-a812-6f6cd89a91c7
_plot(src.twindragon(), p_max)

# ╔═╡ fdc89186-6d96-488b-8bf2-cf97ae18a54a
_plot(src.cantor_dust_non_sym(), p_max)

# ╔═╡ 80907e2f-a22a-4a6d-9fb6-221cecc6fca5
_plot(src.barnsley_fern(), p_max)

# ╔═╡ 5fd68c56-a139-4fb3-952a-a42b05f15f33
_plot(src.cantor_dust(1 / 3, [-1.0, 1.0], 3), p_max)

# ╔═╡ Cell order:
# ╠═49a64b38-6888-11ef-05ba-23348a353dd7
# ╠═1005ccc7-156f-4ea1-a60a-46497e8f144c
# ╠═50acd7c4-0770-4b11-aea4-f852236e4f31
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
# ╠═912bc5b6-a174-4ed0-9fa1-5f4e8de4d5fb
# ╠═7436e59c-1b01-46cf-b6c6-0add22c4b7d4
# ╠═e77088ae-f27b-4b6c-bd1c-881068a95e3b
