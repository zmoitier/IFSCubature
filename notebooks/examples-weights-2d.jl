### A Pluto.jl notebook ###
# v0.20.8

using Markdown
using InteractiveUtils

# ╔═╡ 7970cc66-6bc9-11ef-2104-fb85cdfa2b04
begin
    using Pkg: Pkg
    Pkg.activate(Base.current_project())
    Pkg.instantiate()

    using StaticArrays, GLMakie

    import IFSCubature as src
end

# ╔═╡ dcd6cea7-97fd-4d94-beb6-bc1dcd96b466
begin
    const MAXITER = 2048
    const FONTSIZE = 20
    const SAVEFIG = false

    "Global parameters"
end

# ╔═╡ f3a4b02d-26e8-4824-8ecc-8a3ac7075efb
function _comp_weights(
    sas::src.SelfAffineSet{D,T,N}, pts_cbt_type::String, pts_cbt_max::Int
) where {D,T,N}
    nb_pts = ones(Int, pts_cbt_max)
    weights = ones(Float64, pts_cbt_max)
    for M in 2:pts_cbt_max
        cbt = src.compute_cubature(sas, pts_cbt_type, M; maxiter=MAXITER)
        nb_pts[M] = length(cbt)
        weights[M] = sum(abs.(cbt.weights))
    end

    return (nb_pts, weights)
end

# ╔═╡ 55bc8777-4eae-4b63-829c-3aa24bb719e7
function plot_sum(
    vec_sas::Vector{src.SelfAffineSet{D,T,N}}; pts_cbt_type::String, pts_cbt_max::Int
) where {D,T,N}
    fig = Figure(; size=(600, 400), fontsize=FONTSIZE)
    ax = Axis(fig[1, 1]; yscale=log10)

    for sas in vec_sas
        nb_pts, sum_abs = _comp_weights(sas, pts_cbt_type, pts_cbt_max)
        scatterlines!(ax, nb_pts, sum_abs; linestyle=:dash)
    end

    return fig
end

# ╔═╡ 24535836-fded-4006-9997-b1c90f7968aa
plot_sum(
    [
        src.cantor_dust(1 / 2, [-1.0, 1.0], 2),
        src.cantor_dust(1 / 3, [-1.0, 1.0], 2),
        src.cantor_dust(1 / 4, [-1.0, 1.0], 2),
    ];
    pts_cbt_type="Chebyshev-1",
    pts_cbt_max=16,
)

# ╔═╡ 05dd3e7f-d975-4a2e-a81f-766b29f35911
plot_sum(
    [
        src.cantor_dust(1 / 2, [-1.0, 1.0], 2),
        src.cantor_dust(1 / 3, [-1.0, 1.0], 2),
        src.cantor_dust(1 / 4, [-1.0, 1.0], 2),
    ];
    pts_cbt_type="Equispaced-1",
    pts_cbt_max=16,
)

# ╔═╡ 0abc12ce-d6a4-492a-94da-1e4b20f512ba
plot_sum(
    [src.sierpinski_triangle(), src.sierpinski_triangle(1 / 3)];
    pts_cbt_type="Chebyshev-1",
    pts_cbt_max=16,
)

# ╔═╡ e6d9a47e-33e4-42ee-aabf-019a407843bf
plot_sum(
    [src.sierpinski_triangle(), src.sierpinski_triangle(1 / 3)];
    pts_cbt_type="Equispaced-1",
    pts_cbt_max=16,
)

# ╔═╡ de42d1b9-acab-40c5-be2e-36fa57d35912
plot_sum(
    [src.sierpinski_triangle_fat(2), src.sierpinski_triangle_fat(3)];
    pts_cbt_type="Chebyshev-1",
    pts_cbt_max=16,
)

# ╔═╡ 0b69aec3-3d98-43cf-a6ad-884754bdbf96
plot_sum(
    [src.sierpinski_triangle_fat(2), src.sierpinski_triangle_fat(3)];
    pts_cbt_type="Equispaced-1",
    pts_cbt_max=16,
)

# ╔═╡ 45d43dc7-eb04-4a3a-9a20-004734765c4b
plot_sum(
    [src.vicsek_2d(1 / 3), src.vicsek_2d(1 / 3, 0.4), src.vicsek_2d(1 / 3, π / 4)];
    pts_cbt_type="Chebyshev-1",
    pts_cbt_max=16,
)

# ╔═╡ 41b14410-8956-4637-aade-11ca214a1d17
plot_sum(
    [src.vicsek_2d(1 / 3), src.vicsek_2d(1 / 3, 0.4), src.vicsek_2d(1 / 3, π / 4)];
    pts_cbt_type="Equispaced-1",
    pts_cbt_max=16,
)

# ╔═╡ b8bd21bf-256e-42f2-b9d3-781f662040ad
plot_sum([src.sierpinski_carpet()]; pts_cbt_type="Chebyshev-1", pts_cbt_max=16)

# ╔═╡ f4820474-a6a0-487f-a2c5-096e1a1c0e70
plot_sum([src.sierpinski_carpet()]; pts_cbt_type="Equispaced-1", pts_cbt_max=16)

# ╔═╡ e436f2e5-95c8-4e39-88ec-7fce79d55877
plot_sum([src.koch_snowflake()]; pts_cbt_type="Chebyshev-1", pts_cbt_max=16)

# ╔═╡ 8568010f-d73c-4b4d-99f7-d331fe7cf76f
plot_sum([src.koch_snowflake()]; pts_cbt_type="Equispaced-1", pts_cbt_max=16)

# ╔═╡ cfc8f02d-fca2-44da-be98-461b10152ec7
plot_sum([src.gosper_flowsnake()]; pts_cbt_type="Chebyshev-1", pts_cbt_max=16)

# ╔═╡ 8810e074-f776-4de7-ba2a-7e9438da9e55
plot_sum([src.gosper_flowsnake()]; pts_cbt_type="Equispaced-1", pts_cbt_max=16)

# ╔═╡ 2fc2dd0d-7960-4e7d-906e-dbf2c8bc885a
plot_sum([src.durer_pentagon()]; pts_cbt_type="Chebyshev-1", pts_cbt_max=16)

# ╔═╡ 4347ed03-d6ec-4cf5-9255-2fd6fa61c97d
plot_sum([src.durer_pentagon()]; pts_cbt_type="Equispaced-1", pts_cbt_max=16)

# ╔═╡ 1fd8c1f9-9f31-4fff-a4e9-3d57e26cb0f4
plot_sum([src.fudgeflake()]; pts_cbt_type="Chebyshev-1", pts_cbt_max=16)

# ╔═╡ 5e13c97f-abfa-45e1-af43-b16deb16565d
plot_sum([src.fudgeflake()]; pts_cbt_type="Equispaced-1", pts_cbt_max=16)

# ╔═╡ 7c843437-d0eb-44d7-8ffd-ab073ed95a57
plot_sum([src.heighway_dragon()]; pts_cbt_type="Chebyshev-1", pts_cbt_max=16)

# ╔═╡ c3400ca4-5cae-4530-8390-5a80979ac80d
plot_sum([src.levy_dragon()]; pts_cbt_type="Chebyshev-1", pts_cbt_max=16)

# ╔═╡ c4822134-b393-4986-afcd-4a42c5216523
plot_sum([src.levy_dragon()]; pts_cbt_type="Equispaced-1", pts_cbt_max=16)

# ╔═╡ b7e769bc-4f97-46fc-a755-bea6707545d6
plot_sum([src.terdragon()]; pts_cbt_type="Chebyshev-1", pts_cbt_max=16)

# ╔═╡ aff2d174-bca3-4fad-af21-7cfa766f3796
plot_sum([src.terdragon()]; pts_cbt_type="Equispaced-1", pts_cbt_max=16)

# ╔═╡ b4fae1ca-ea2c-4009-b3c0-852b33189418
plot_sum([src.twindragon()]; pts_cbt_type="Chebyshev-1", pts_cbt_max=16)

# ╔═╡ cbc2786b-6200-4e7f-91a7-ee0fab5ad47b
plot_sum([src.twindragon()]; pts_cbt_type="Equispaced-1", pts_cbt_max=16)

# ╔═╡ 6411eb71-9f7c-48bd-9c2f-b2cab7ab7184
plot_sum([src.brick_2d()]; pts_cbt_type="Chebyshev-1", pts_cbt_max=16)

# ╔═╡ 197f1284-25d1-421d-8efd-ac4d13e503df
plot_sum([src.brick_2d()]; pts_cbt_type="Equispaced-1", pts_cbt_max=16)

# ╔═╡ 2314de5c-db07-4d10-b17e-a8cba8d798e1
plot_sum([src.cantor_dust_non_sym()]; pts_cbt_type="Chebyshev-1", pts_cbt_max=16)

# ╔═╡ 50799326-425e-4fca-88e2-8c78384aaf7a
plot_sum([src.cantor_dust_non_sym()]; pts_cbt_type="Equispaced-1", pts_cbt_max=16)

# ╔═╡ acc59a6f-a5f9-4ac3-9b26-ba5a2be5f780
plot_sum([src.barnsley_fern()]; pts_cbt_type="Chebyshev-1", pts_cbt_max=16)

# ╔═╡ 67447ed3-9cb4-4507-810c-f12bf6d4eb64
plot_sum([src.barnsley_fern()]; pts_cbt_type="Equispaced-1", pts_cbt_max=16)

# ╔═╡ 7b5d13fe-8911-48c8-8588-b43dbab5963a
function plot_chaos_game!(ax, sas::src.SelfAffineSet{2,Float64,4}, nb_pts::Int, α::Real)
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

    v = collect(0:(n - 1)) .* h
    heatmap!(
        ax,
        _min[1] .+ v,
        _min[2] .+ v,
        f;
        colormap=Reverse(:grays),
        colorrange=(0, 1),
        alpha=α,
    )

    return nothing
end

# ╔═╡ 5cb0a141-15bd-4dbd-a833-a1415730db2f
function _get_limits(box::src.HyperBox{2,T,4}) where {T}
    vertices = src.vertices(box)
    xmin, xmax = extrema(v[1] for v in vertices)
    ymin, ymax = extrema(v[2] for v in vertices)
    h, w = 0.05(xmax - xmin), 0.05 * (ymax - ymin)
    return (xmin - h, xmax + h, ymin - h, ymax + h)
end

# ╔═╡ cfcd7ec2-6c6b-4400-b9fd-0bc8df2059ac
function plot_weights_2d(;
    sas::src.SelfAffineSet{2,T,4},
    pts_chaos_nb::Int,
    α_attractor::Real,
    pts_cbt_type::String,
    pts_cbt_nb::Int,
    α_weights::Real,
) where {T}
    fig = Figure(; size=(600, 600), fontsize=FONTSIZE)
    ax = Axis(fig[1, 1]; aspect=1, xlabel=L"x", ylabel=L"y")

    plot_chaos_game!(ax, sas, pts_chaos_nb, α_attractor)

    cbt = src.compute_cubature(sas, pts_cbt_type, pts_cbt_nb; maxiter=MAXITER)
    idx_pos = findall(>(0), cbt.weights)
    idx_neg = findall(<(0), cbt.weights)

    colors = Makie.to_colormap(:tab10)
    scatter!(ax, cbt.points[idx_pos]; color=(colors[3], α_weights))
    scatter!(ax, cbt.points[idx_neg]; color=(colors[4], α_weights))

    xmin, xmax, ymin, ymax = _get_limits(sas.bounding_box)
    limits!(ax, xmin, xmax, ymin, ymax)

    if SAVEFIG
        save("$name-weights-2d.pdf", fig)
    end

    return fig
end

# ╔═╡ 0420f787-d417-42d7-8a4e-c4b5f7078717
plot_weights_2d(;
    sas=src.cantor_dust(1 / 3, [-1.0, 1.0], 2),
    pts_chaos_nb=200_000,
    α_attractor=0.25,
    pts_cbt_type="Chebyshev-1",
    pts_cbt_nb=16,
    α_weights=0.75,
)

# ╔═╡ 6136f605-38fd-4046-95f3-a16400dc4d34
plot_weights_2d(;
    sas=src.sierpinski_triangle(),
    pts_chaos_nb=200_000,
    α_attractor=0.25,
    pts_cbt_type="Chebyshev-1",
    pts_cbt_nb=16,
    α_weights=0.75,
)

# ╔═╡ 3fbf1a9c-5554-4594-b260-98abd01f56db
plot_weights_2d(;
    sas=src.sierpinski_triangle_fat(2),
    pts_chaos_nb=200_000,
    α_attractor=0.25,
    pts_cbt_type="Chebyshev-1",
    pts_cbt_nb=16,
    α_weights=0.75,
)

# ╔═╡ 22fd1def-fb18-4910-b1c4-16fd54003e8a
plot_weights_2d(;
    sas=src.vicsek_2d(1 / 3),
    pts_chaos_nb=200_000,
    α_attractor=0.25,
    pts_cbt_type="Chebyshev-1",
    pts_cbt_nb=16,
    α_weights=0.75,
)

# ╔═╡ ee205c3f-5514-42a5-8eeb-a4e1447229ad
plot_weights_2d(;
    sas=src.vicsek_2d(1 / 3, 0.4),
    pts_chaos_nb=200_000,
    α_attractor=0.25,
    pts_cbt_type="Chebyshev-1",
    pts_cbt_nb=16,
    α_weights=0.75,
)

# ╔═╡ c23fea07-1737-44c9-8e0f-4d4dfbae6cfd
plot_weights_2d(;
    sas=src.vicsek_2d(1 / 3, π / 4),
    pts_chaos_nb=200_000,
    α_attractor=0.25,
    pts_cbt_type="Chebyshev-1",
    pts_cbt_nb=16,
    α_weights=0.75,
)

# ╔═╡ bb4ee489-81aa-4ae8-af03-2685cec5219a
plot_weights_2d(;
    sas=src.sierpinski_carpet(),
    pts_chaos_nb=200_000,
    α_attractor=0.25,
    pts_cbt_type="Chebyshev-1",
    pts_cbt_nb=16,
    α_weights=0.75,
)

# ╔═╡ cc17ee44-9769-4c91-b08c-af138a4d29b7
plot_weights_2d(;
    sas=src.koch_snowflake(),
    pts_chaos_nb=200_000,
    α_attractor=0.25,
    pts_cbt_type="Chebyshev-1",
    pts_cbt_nb=16,
    α_weights=0.75,
)

# ╔═╡ c8d3184c-cac1-46e2-9f97-924aae266444
plot_weights_2d(;
    sas=src.gosper_flowsnake(),
    pts_chaos_nb=200_000,
    α_attractor=0.25,
    pts_cbt_type="Chebyshev-1",
    pts_cbt_nb=16,
    α_weights=0.75,
)

# ╔═╡ a35bb1c9-d044-46bb-806a-4434013881b1
plot_weights_2d(;
    sas=src.durer_pentagon(),
    pts_chaos_nb=200_000,
    α_attractor=0.25,
    pts_cbt_type="Chebyshev-1",
    pts_cbt_nb=16,
    α_weights=0.75,
)

# ╔═╡ 565d6e26-3926-46c4-a531-cf58fcd48d56
plot_weights_2d(;
    sas=src.fudgeflake(),
    pts_chaos_nb=200_000,
    α_attractor=0.25,
    pts_cbt_type="Chebyshev-1",
    pts_cbt_nb=16,
    α_weights=0.75,
)

# ╔═╡ de33fa02-2141-4c86-8eac-56624c2e26d3
plot_weights_2d(;
    sas=src.heighway_dragon(),
    pts_chaos_nb=200_000,
    α_attractor=0.25,
    pts_cbt_type="Chebyshev-1",
    pts_cbt_nb=16,
    α_weights=0.75,
)

# ╔═╡ 1e60655e-7c4f-4639-9d3b-8509056bfa1c
plot_weights_2d(;
    sas=src.levy_dragon(),
    pts_chaos_nb=200_000,
    α_attractor=0.25,
    pts_cbt_type="Chebyshev-1",
    pts_cbt_nb=16,
    α_weights=0.75,
)

# ╔═╡ 903b2471-753d-4af2-a2b4-8b3b56ad8660
plot_weights_2d(;
    sas=src.terdragon(),
    pts_chaos_nb=200_000,
    α_attractor=0.75,
    pts_cbt_type="Chebyshev-1",
    pts_cbt_nb=16,
    α_weights=0.75,
)

# ╔═╡ e3eb8f48-c922-43f3-9b8f-37ff282fdaa3
plot_weights_2d(;
    sas=src.twindragon(),
    pts_chaos_nb=200_000,
    α_attractor=0.75,
    pts_cbt_type="Chebyshev-1",
    pts_cbt_nb=16,
    α_weights=0.75,
)

# ╔═╡ d06989a2-d646-4019-a4d2-16e242c1b0ed
plot_weights_2d(;
    sas=src.brick_2d(),
    pts_chaos_nb=200_000,
    α_attractor=0.75,
    pts_cbt_type="Chebyshev-1",
    pts_cbt_nb=16,
    α_weights=0.75,
)

# ╔═╡ 14d64bca-e80a-4214-be64-08a6ae078241
plot_weights_2d(;
    sas=src.cantor_dust_non_sym(),
    pts_chaos_nb=200_000,
    α_attractor=0.75,
    pts_cbt_type="Chebyshev-1",
    pts_cbt_nb=16,
    α_weights=0.75,
)

# ╔═╡ da548fdb-98eb-4f76-8d94-d7063d7d7f9c
plot_weights_2d(;
    sas=src.barnsley_fern(),
    pts_chaos_nb=200_000,
    α_attractor=0.75,
    pts_cbt_type="Chebyshev-1",
    pts_cbt_nb=16,
    α_weights=0.75,
)

# ╔═╡ c20338c7-122e-43a1-b937-2da72f008a02
function plot_weights_3d(;
    sas::src.SelfAffineSet{2,T,4},
    pts_chaos_nb::Int,
    α_attractor::Real,
    pts_cbt_type::String,
    pts_cbt_nb::Int,
    α_weights::Real,
) where {T}
    fig = Figure(; size=(800, 600), fontsize=FONTSIZE)
    ax = Axis3(fig[1, 1]; xlabel=L"x", ylabel=L"y", zlabel=L"|w|")

    plot_chaos_game!(ax, sas, pts_chaos_nb, α_attractor)

    cbt = src.compute_cubature(sas, pts_cbt_type, pts_cbt_nb; maxiter=MAXITER)
    idx_pos = findall(>(0), cbt.weights)
    idx_neg = findall(<(0), cbt.weights)

    colors = Makie.to_colormap(:tab10)
    for (idx, ic) in zip((findall(>(0), cbt.weights), findall(<(0), cbt.weights)), (3, 4))
        for (x, w) in zip(cbt.points[idx], cbt.weights[idx])
            lines!(
                ax,
                [x[1], x[1]],
                [x[2], x[2]],
                [0.0, abs.(w)];
                color=(colors[ic], 0.75),
                linestyle=:dash,
            )
            scatter!(ax, x[1], x[2], abs.(w); color=(colors[ic], α_weights))
        end
    end

    xmin, xmax, ymin, ymax = _get_limits(sas.bounding_box)
    limits!(ax, xmin, xmax, ymin, ymax, 0.0, maximum(abs.(cbt.weights)))

    if SAVEFIG
        save("$name-weights-3d.pdf", fig)
    end

    return fig
end

# ╔═╡ ec30462f-d3dc-418f-962d-c109fc9f50f1
plot_weights_3d(;
    sas=src.cantor_dust(1 / 3, [-1.0, 1.0], 2),
    pts_chaos_nb=200_000,
    α_attractor=0.75,
    pts_cbt_type="Chebyshev-1",
    pts_cbt_nb=16,
    α_weights=0.75,
)

# ╔═╡ 83ea3c34-48b9-4eba-8982-af7cb470df8c
plot_weights_3d(;
    sas=src.sierpinski_triangle(),
    pts_chaos_nb=200_000,
    α_attractor=0.75,
    pts_cbt_type="Chebyshev-1",
    pts_cbt_nb=16,
    α_weights=0.75,
)

# ╔═╡ 911f7450-5e3d-4b2a-9d74-ac128292e1aa
plot_weights_3d(;
    sas=src.sierpinski_triangle_fat(2),
    pts_chaos_nb=200_000,
    α_attractor=0.75,
    pts_cbt_type="Chebyshev-1",
    pts_cbt_nb=16,
    α_weights=0.75,
)

# ╔═╡ 218e71d7-93dd-4b8b-bb5a-ccf39555bbf6
plot_weights_3d(;
    sas=src.vicsek_2d(1 / 3),
    pts_chaos_nb=200_000,
    α_attractor=0.75,
    pts_cbt_type="Chebyshev-1",
    pts_cbt_nb=16,
    α_weights=0.75,
)

# ╔═╡ e38481e0-c432-4220-beaa-f7912e0c56f1
plot_weights_3d(;
    sas=src.vicsek_2d(1 / 3, 0.4),
    pts_chaos_nb=200_000,
    α_attractor=0.75,
    pts_cbt_type="Chebyshev-1",
    pts_cbt_nb=16,
    α_weights=0.75,
)

# ╔═╡ f2c95728-c94f-48d3-a7a8-a8816a8e0cca
plot_weights_3d(;
    sas=src.vicsek_2d(1 / 3, π / 4),
    pts_chaos_nb=200_000,
    α_attractor=0.75,
    pts_cbt_type="Chebyshev-1",
    pts_cbt_nb=16,
    α_weights=0.75,
)

# ╔═╡ 4ab0c8f4-b370-4107-9c5b-98bc039c3df9
plot_weights_3d(;
    sas=src.sierpinski_carpet(),
    pts_chaos_nb=200_000,
    α_attractor=0.75,
    pts_cbt_type="Chebyshev-1",
    pts_cbt_nb=16,
    α_weights=0.75,
)

# ╔═╡ 623f596a-3491-4483-982d-163eade449a3
plot_weights_3d(;
    sas=src.koch_snowflake(),
    pts_chaos_nb=200_000,
    α_attractor=0.75,
    pts_cbt_type="Chebyshev-1",
    pts_cbt_nb=16,
    α_weights=0.75,
)

# ╔═╡ 44af8667-d453-4db1-a305-f63629e1281c
plot_weights_3d(;
    sas=src.gosper_flowsnake(),
    pts_chaos_nb=200_000,
    α_attractor=0.75,
    pts_cbt_type="Chebyshev-1",
    pts_cbt_nb=16,
    α_weights=0.75,
)

# ╔═╡ 14d645d9-a1f1-4a35-8fa1-27965f9fab76
plot_weights_3d(;
    sas=src.durer_pentagon(),
    pts_chaos_nb=100_000,
    α_attractor=0.75,
    pts_cbt_type="Chebyshev-1",
    pts_cbt_nb=16,
    α_weights=0.75,
)

# ╔═╡ ae142920-e802-4842-8d87-2fcd81b8317d
plot_weights_3d(;
    sas=src.fudgeflake(),
    pts_chaos_nb=100_000,
    α_attractor=0.75,
    pts_cbt_type="Chebyshev-1",
    pts_cbt_nb=16,
    α_weights=0.75,
)

# ╔═╡ a4c7afa7-75c2-431c-a7ff-c874d17248d3
plot_weights_3d(;
    sas=src.heighway_dragon(),
    pts_chaos_nb=200_000,
    α_attractor=0.75,
    pts_cbt_type="Chebyshev-1",
    pts_cbt_nb=16,
    α_weights=0.75,
)

# ╔═╡ 3a1f37b7-a8ec-4371-a37f-468587168240
plot_weights_3d(;
    sas=src.levy_dragon(),
    pts_chaos_nb=200_000,
    α_attractor=0.75,
    pts_cbt_type="Chebyshev-1",
    pts_cbt_nb=16,
    α_weights=0.75,
)

# ╔═╡ 002ad1a3-f7e3-42a7-8aed-aa8b40294fe2
plot_weights_3d(;
    sas=src.terdragon(),
    pts_chaos_nb=200_000,
    α_attractor=0.75,
    pts_cbt_type="Chebyshev-1",
    pts_cbt_nb=16,
    α_weights=0.75,
)

# ╔═╡ 2ffee5ea-5191-4759-9785-72d37058e641
plot_weights_3d(;
    sas=src.twindragon(),
    pts_chaos_nb=200_000,
    α_attractor=0.75,
    pts_cbt_type="Chebyshev-1",
    pts_cbt_nb=16,
    α_weights=0.75,
)

# ╔═╡ bcaed2e5-57ca-4ade-98de-0a6e30f2d780
plot_weights_3d(;
    sas=src.brick_2d(),
    pts_chaos_nb=200_000,
    α_attractor=0.75,
    pts_cbt_type="Chebyshev-1",
    pts_cbt_nb=16,
    α_weights=0.75,
)

# ╔═╡ 34ac624f-8c8d-4d70-9f6f-33cafc3b20cd
plot_weights_3d(;
    sas=src.cantor_dust_non_sym(),
    pts_chaos_nb=200_000,
    α_attractor=0.75,
    pts_cbt_type="Chebyshev-1",
    pts_cbt_nb=16,
    α_weights=0.75,
)

# ╔═╡ dae1ff12-5441-43b4-8df4-0d60dd6732bc
plot_weights_3d(;
    sas=src.barnsley_fern(),
    pts_chaos_nb=200_000,
    α_attractor=0.75,
    pts_cbt_type="Chebyshev-1",
    pts_cbt_nb=16,
    α_weights=0.75,
)

# ╔═╡ Cell order:
# ╠═7970cc66-6bc9-11ef-2104-fb85cdfa2b04
# ╠═dcd6cea7-97fd-4d94-beb6-bc1dcd96b466
# ╠═24535836-fded-4006-9997-b1c90f7968aa
# ╠═05dd3e7f-d975-4a2e-a81f-766b29f35911
# ╠═0420f787-d417-42d7-8a4e-c4b5f7078717
# ╠═ec30462f-d3dc-418f-962d-c109fc9f50f1
# ╠═0abc12ce-d6a4-492a-94da-1e4b20f512ba
# ╠═e6d9a47e-33e4-42ee-aabf-019a407843bf
# ╠═6136f605-38fd-4046-95f3-a16400dc4d34
# ╠═83ea3c34-48b9-4eba-8982-af7cb470df8c
# ╠═de42d1b9-acab-40c5-be2e-36fa57d35912
# ╠═0b69aec3-3d98-43cf-a6ad-884754bdbf96
# ╠═3fbf1a9c-5554-4594-b260-98abd01f56db
# ╠═911f7450-5e3d-4b2a-9d74-ac128292e1aa
# ╠═45d43dc7-eb04-4a3a-9a20-004734765c4b
# ╠═41b14410-8956-4637-aade-11ca214a1d17
# ╠═22fd1def-fb18-4910-b1c4-16fd54003e8a
# ╠═218e71d7-93dd-4b8b-bb5a-ccf39555bbf6
# ╠═ee205c3f-5514-42a5-8eeb-a4e1447229ad
# ╠═e38481e0-c432-4220-beaa-f7912e0c56f1
# ╠═c23fea07-1737-44c9-8e0f-4d4dfbae6cfd
# ╠═f2c95728-c94f-48d3-a7a8-a8816a8e0cca
# ╠═b8bd21bf-256e-42f2-b9d3-781f662040ad
# ╠═f4820474-a6a0-487f-a2c5-096e1a1c0e70
# ╠═bb4ee489-81aa-4ae8-af03-2685cec5219a
# ╠═4ab0c8f4-b370-4107-9c5b-98bc039c3df9
# ╠═e436f2e5-95c8-4e39-88ec-7fce79d55877
# ╠═8568010f-d73c-4b4d-99f7-d331fe7cf76f
# ╠═cc17ee44-9769-4c91-b08c-af138a4d29b7
# ╠═623f596a-3491-4483-982d-163eade449a3
# ╠═cfc8f02d-fca2-44da-be98-461b10152ec7
# ╠═8810e074-f776-4de7-ba2a-7e9438da9e55
# ╠═c8d3184c-cac1-46e2-9f97-924aae266444
# ╠═44af8667-d453-4db1-a305-f63629e1281c
# ╠═2fc2dd0d-7960-4e7d-906e-dbf2c8bc885a
# ╠═4347ed03-d6ec-4cf5-9255-2fd6fa61c97d
# ╠═a35bb1c9-d044-46bb-806a-4434013881b1
# ╠═14d645d9-a1f1-4a35-8fa1-27965f9fab76
# ╠═1fd8c1f9-9f31-4fff-a4e9-3d57e26cb0f4
# ╠═5e13c97f-abfa-45e1-af43-b16deb16565d
# ╠═565d6e26-3926-46c4-a531-cf58fcd48d56
# ╠═ae142920-e802-4842-8d87-2fcd81b8317d
# ╠═7c843437-d0eb-44d7-8ffd-ab073ed95a57
# ╠═de33fa02-2141-4c86-8eac-56624c2e26d3
# ╠═a4c7afa7-75c2-431c-a7ff-c874d17248d3
# ╠═c3400ca4-5cae-4530-8390-5a80979ac80d
# ╠═c4822134-b393-4986-afcd-4a42c5216523
# ╠═1e60655e-7c4f-4639-9d3b-8509056bfa1c
# ╠═3a1f37b7-a8ec-4371-a37f-468587168240
# ╠═b7e769bc-4f97-46fc-a755-bea6707545d6
# ╠═aff2d174-bca3-4fad-af21-7cfa766f3796
# ╠═903b2471-753d-4af2-a2b4-8b3b56ad8660
# ╠═002ad1a3-f7e3-42a7-8aed-aa8b40294fe2
# ╠═b4fae1ca-ea2c-4009-b3c0-852b33189418
# ╠═cbc2786b-6200-4e7f-91a7-ee0fab5ad47b
# ╠═e3eb8f48-c922-43f3-9b8f-37ff282fdaa3
# ╠═2ffee5ea-5191-4759-9785-72d37058e641
# ╠═6411eb71-9f7c-48bd-9c2f-b2cab7ab7184
# ╠═197f1284-25d1-421d-8efd-ac4d13e503df
# ╠═d06989a2-d646-4019-a4d2-16e242c1b0ed
# ╠═bcaed2e5-57ca-4ade-98de-0a6e30f2d780
# ╠═2314de5c-db07-4d10-b17e-a8cba8d798e1
# ╠═50799326-425e-4fca-88e2-8c78384aaf7a
# ╠═14d64bca-e80a-4214-be64-08a6ae078241
# ╠═34ac624f-8c8d-4d70-9f6f-33cafc3b20cd
# ╠═acc59a6f-a5f9-4ac3-9b26-ba5a2be5f780
# ╠═67447ed3-9cb4-4507-810c-f12bf6d4eb64
# ╠═da548fdb-98eb-4f76-8d94-d7063d7d7f9c
# ╠═dae1ff12-5441-43b4-8df4-0d60dd6732bc
# ╠═55bc8777-4eae-4b63-829c-3aa24bb719e7
# ╠═f3a4b02d-26e8-4824-8ecc-8a3ac7075efb
# ╠═cfcd7ec2-6c6b-4400-b9fd-0bc8df2059ac
# ╠═c20338c7-122e-43a1-b937-2da72f008a02
# ╠═7b5d13fe-8911-48c8-8588-b43dbab5963a
# ╠═5cb0a141-15bd-4dbd-a833-a1415730db2f
