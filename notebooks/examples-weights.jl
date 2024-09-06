### A Pluto.jl notebook ###
# v0.19.46

using Markdown
using InteractiveUtils

# ╔═╡ 7970cc66-6bc9-11ef-2104-fb85cdfa2b04
begin
    using Pkg: Pkg
    Pkg.activate(Base.current_project())
    Pkg.instantiate()

    using CairoMakie

    import IFSCubature as src
end

# ╔═╡ dcd6cea7-97fd-4d94-beb6-bc1dcd96b466
begin
    #! Type of points
    # const POINTTYPE = "equispaced"
    const POINTTYPE = "Chebyshev-1"
    # const POINTTYPE = "Chebyshev-2"
    # const POINTTYPE = "Gauss-Legendre"
    # const POINTTYPE = "Gauss-Lobatto"

    const MAXITER = 2048

    const FONTSIZE = 20

    "Global parameters"
end

# ╔═╡ f6ab7b37-ae52-4998-a8f6-295fc67360cd
function _save_weights(sas::src.SelfAffineSet{D,T,N}, nb_pts_max::Int) where {D,T,N}
    nb_pts = ones(Int, nb_pts_max)
    sum_abs = ones(Float64, nb_pts_max)
    for M in 2:nb_pts_max
        cbt = src.compute_cubature(sas, POINTTYPE, M; maxiter=MAXITER)
        nb_pts[M] = length(cbt)
        sum_abs[M] = sum(abs.(cbt.weights))
    end

    return (nb_pts, sum_abs)
end

# ╔═╡ 4c984422-dbaa-4414-98d8-74fdcbbfd954
function clean_name(name::String)
    return replace(name, "-" => " ")
end

# ╔═╡ 238dbcfd-bc15-443b-935c-0960fe1d2b0b
function plot_weights_1d(;
    sas::src.SelfAffineSet{1,T,1}, nb_pts_cbt::Int, nb_refine::Int
) where {T}
    cbt = src.compute_cubature(sas, POINTTYPE, nb_pts_cbt; maxiter=MAXITER)
    w_lim = extrema(cbt.weights)
    W = maximum(abs.(w_lim))

    dx, dy = 0.025, 0.05 * (w_lim[2] - w_lim[1])
    X = (-dx, 1 + dx)
    Y = (w_lim[1] - dy, w_lim[2] + dy)

    boxes = [sas.bounding_box]
    for _ in 1:nb_refine
        boxes = [S(box) for S in sas.ifs for box in boxes]
    end

    fig = Figure(; fontsize=FONTSIZE)
    ax = Axis(
        fig[1, 1];
        title=L"%$(clean_name(sas.name)) and $ M = %$nb_pts_cbt $",
        xlabel=L"x",
        ylabel=L"w",
    )

    limits!(ax, X, Y)

    for box in boxes
        a, b = src.vertices(box)
        poly!(ax, [(a, Y[1]), (b, Y[1]), (b, Y[2]), (a, Y[2])]; color=(:black, 0.3))
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

    return fig
end

# ╔═╡ a3bbe18f-c413-4435-9d8f-6e559615fe26
plot_weights_1d(; sas=src.cantor_set(1 / 3, [0.0, 1.0]), nb_pts_cbt=127, nb_refine=5)

# ╔═╡ 17239ef0-6170-4f11-99dd-d7cd63b49775
function plot_weights_2d(;
    sas::src.SelfAffineSet{2,T,4}, nb_pts_cbt::Int, nb_refine::Int
) where {T}
    cbt = src.compute_cubature(sas, POINTTYPE, nb_pts_cbt; maxiter=MAXITER)
    W = maximum(abs.(cbt.weights))

    pf = [src.vertices(sas.bounding_box)]
    for _ in 1:nb_refine
        pf = [S.(part) for S in sas.ifs for part in pf]
    end

    fig = Figure(; fontsize=FONTSIZE)

    M2 = nb_pts_cbt * nb_pts_cbt
    ax = Axis(
        fig[1, 1];
        title=L"%$(clean_name(sas.name)) and $ M = %$M2 $",
        xlabel=L"x",
        ylabel=L"y",
    )

    dx = 0.025
    # limits!(ax, (-1 - dx, 1 + dx), (-1 - dx, 1 + dx))

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

    Colorbar(fig[1, 2], sc; label=L"w")

    return fig
end

# ╔═╡ 10adbc06-3ccc-42a9-a977-881fb142a6f8
plot_weights_2d(; sas=src.cantor_dust(1 / 3, [-1.0, 1.0], 2), nb_pts_cbt=10, nb_refine=3)

# ╔═╡ 6136f605-38fd-4046-95f3-a16400dc4d34
plot_weights_2d(; sas=src.sierpinski_triangle(), nb_pts_cbt=10, nb_refine=3)

# ╔═╡ 3fbf1a9c-5554-4594-b260-98abd01f56db
plot_weights_2d(; sas=src.sierpinski_triangle_fat(2), nb_pts_cbt=10, nb_refine=3)

# ╔═╡ 22fd1def-fb18-4910-b1c4-16fd54003e8a
plot_weights_2d(; sas=src.vicsek_2d(1 / 3), nb_pts_cbt=10, nb_refine=3)

# ╔═╡ 218e71d7-93dd-4b8b-bb5a-ccf39555bbf6
plot_weights_2d(; sas=src.vicsek_2d(1 / 3, π / 4), nb_pts_cbt=10, nb_refine=3)

# ╔═╡ bb4ee489-81aa-4ae8-af03-2685cec5219a
plot_weights_2d(; sas=src.sierpinski_carpet(), nb_pts_cbt=10, nb_refine=3)

# ╔═╡ cc17ee44-9769-4c91-b08c-af138a4d29b7
plot_weights_2d(; sas=src.koch_snowflake(), nb_pts_cbt=10, nb_refine=3)

# ╔═╡ c8d3184c-cac1-46e2-9f97-924aae266444
plot_weights_2d(; sas=src.gosper_flowsnake(), nb_pts_cbt=10, nb_refine=3)

# ╔═╡ d06989a2-d646-4019-a4d2-16e242c1b0ed
plot_weights_2d(; sas=src.brick_2d(), nb_pts_cbt=10, nb_refine=3)

# ╔═╡ 14d64bca-e80a-4214-be64-08a6ae078241
plot_weights_2d(; sas=src.cantor_dust_non_sym(), nb_pts_cbt=10, nb_refine=3)

# ╔═╡ da548fdb-98eb-4f76-8d94-d7063d7d7f9c
plot_weights_2d(; sas=src.barnsley_fern(), nb_pts_cbt=10, nb_refine=3)

# ╔═╡ 55bc8777-4eae-4b63-829c-3aa24bb719e7
function plot_sum(; sas::src.SelfAffineSet{D,T,N}, nb_pts_max::Int) where {D,T,N}
    nb_pts, sum_abs = _save_weights(sas, nb_pts_max)

    fig = Figure(; fontsize=FONTSIZE)
    ax = Axis(
        fig[1, 1];
        title=L"%$(clean_name(sas.name))$$",
        xlabel=L"number of cubature points$$",
        ylabel=L"|\mathbf{w}|_1",
    )

    scatterlines!(ax, nb_pts, sum_abs; linestyle=:dash)

    return fig
end

# ╔═╡ 0f6f889e-6a7f-4d0a-b32e-d84eb6859c8d
plot_sum(; sas=src.cantor_set(1 / 3, [0.0, 1.0]), nb_pts_max=64)

# ╔═╡ 24535836-fded-4006-9997-b1c90f7968aa
plot_sum(; sas=src.cantor_dust(1 / 3, [-1.0, 1.0], 2), nb_pts_max=20)

# ╔═╡ ab254274-e56c-40a4-ab58-4accdbb8b07f
plot_sum(; sas=src.sierpinski_triangle(), nb_pts_max=20)

# ╔═╡ 1dac0461-34b0-4ad9-a3f7-d2532c86f489
plot_sum(; sas=src.sierpinski_triangle_fat(2), nb_pts_max=20)

# ╔═╡ da36f44c-0f74-43f7-a2a0-9dddba121692
plot_sum(; sas=src.vicsek_2d(1 / 3), nb_pts_max=20)

# ╔═╡ 50a34f8e-7726-4b83-8e22-9d989ae829bd
plot_sum(; sas=src.vicsek_2d(1 / 3, π / 4), nb_pts_max=20)

# ╔═╡ abf659ae-6659-4a55-9169-7e8c0eba788f
plot_sum(; sas=src.sierpinski_carpet(), nb_pts_max=20)

# ╔═╡ 209b79d9-3bc4-4a68-857e-327e91cbbe2e
plot_sum(; sas=src.koch_snowflake(), nb_pts_max=20)

# ╔═╡ 44af8667-d453-4db1-a305-f63629e1281c
plot_sum(; sas=src.gosper_flowsnake(), nb_pts_max=20)

# ╔═╡ bcaed2e5-57ca-4ade-98de-0a6e30f2d780
plot_sum(; sas=src.brick_2d(), nb_pts_max=20)

# ╔═╡ 34ac624f-8c8d-4d70-9f6f-33cafc3b20cd
plot_sum(; sas=src.cantor_dust_non_sym(), nb_pts_max=20)

# ╔═╡ dae1ff12-5441-43b4-8df4-0d60dd6732bc
plot_sum(; sas=src.barnsley_fern(), nb_pts_max=20)

# ╔═╡ Cell order:
# ╠═7970cc66-6bc9-11ef-2104-fb85cdfa2b04
# ╠═dcd6cea7-97fd-4d94-beb6-bc1dcd96b466
# ╠═a3bbe18f-c413-4435-9d8f-6e559615fe26
# ╠═0f6f889e-6a7f-4d0a-b32e-d84eb6859c8d
# ╠═10adbc06-3ccc-42a9-a977-881fb142a6f8
# ╠═24535836-fded-4006-9997-b1c90f7968aa
# ╠═6136f605-38fd-4046-95f3-a16400dc4d34
# ╠═ab254274-e56c-40a4-ab58-4accdbb8b07f
# ╠═3fbf1a9c-5554-4594-b260-98abd01f56db
# ╠═1dac0461-34b0-4ad9-a3f7-d2532c86f489
# ╠═22fd1def-fb18-4910-b1c4-16fd54003e8a
# ╠═da36f44c-0f74-43f7-a2a0-9dddba121692
# ╠═218e71d7-93dd-4b8b-bb5a-ccf39555bbf6
# ╠═50a34f8e-7726-4b83-8e22-9d989ae829bd
# ╠═bb4ee489-81aa-4ae8-af03-2685cec5219a
# ╠═abf659ae-6659-4a55-9169-7e8c0eba788f
# ╠═cc17ee44-9769-4c91-b08c-af138a4d29b7
# ╠═209b79d9-3bc4-4a68-857e-327e91cbbe2e
# ╠═c8d3184c-cac1-46e2-9f97-924aae266444
# ╠═44af8667-d453-4db1-a305-f63629e1281c
# ╠═d06989a2-d646-4019-a4d2-16e242c1b0ed
# ╠═bcaed2e5-57ca-4ade-98de-0a6e30f2d780
# ╠═14d64bca-e80a-4214-be64-08a6ae078241
# ╠═34ac624f-8c8d-4d70-9f6f-33cafc3b20cd
# ╠═da548fdb-98eb-4f76-8d94-d7063d7d7f9c
# ╠═dae1ff12-5441-43b4-8df4-0d60dd6732bc
# ╠═238dbcfd-bc15-443b-935c-0960fe1d2b0b
# ╠═17239ef0-6170-4f11-99dd-d7cd63b49775
# ╠═55bc8777-4eae-4b63-829c-3aa24bb719e7
# ╠═f6ab7b37-ae52-4998-a8f6-295fc67360cd
# ╠═4c984422-dbaa-4414-98d8-74fdcbbfd954
