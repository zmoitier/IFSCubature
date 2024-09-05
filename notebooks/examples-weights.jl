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

# ╔═╡ Cell order:
# ╠═7970cc66-6bc9-11ef-2104-fb85cdfa2b04
# ╠═dcd6cea7-97fd-4d94-beb6-bc1dcd96b466
# ╠═a3bbe18f-c413-4435-9d8f-6e559615fe26
# ╠═10adbc06-3ccc-42a9-a977-881fb142a6f8
# ╠═6136f605-38fd-4046-95f3-a16400dc4d34
# ╠═3fbf1a9c-5554-4594-b260-98abd01f56db
# ╠═22fd1def-fb18-4910-b1c4-16fd54003e8a
# ╠═218e71d7-93dd-4b8b-bb5a-ccf39555bbf6
# ╠═bb4ee489-81aa-4ae8-af03-2685cec5219a
# ╠═cc17ee44-9769-4c91-b08c-af138a4d29b7
# ╠═c8d3184c-cac1-46e2-9f97-924aae266444
# ╠═d06989a2-d646-4019-a4d2-16e242c1b0ed
# ╠═14d64bca-e80a-4214-be64-08a6ae078241
# ╠═238dbcfd-bc15-443b-935c-0960fe1d2b0b
# ╠═17239ef0-6170-4f11-99dd-d7cd63b49775
# ╠═4c984422-dbaa-4414-98d8-74fdcbbfd954
