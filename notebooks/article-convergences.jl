### A Pluto.jl notebook ###
# v0.19.46

using Markdown
using InteractiveUtils

# ╔═╡ fe657a18-63c0-11ef-0959-15170cb4d30a
begin
    using Pkg: Pkg
    Pkg.activate(Base.current_project())
    Pkg.instantiate()

    using Printf, LinearAlgebra, StaticArrays, DataStructures, TOML, CairoMakie

    import IFSCubature as src
end

# ╔═╡ b5d45ad9-3996-401f-96e5-78aba17e8686
begin
    #! Type of points
    # const POINTTYPE = "equispaced"
    const POINTTYPE = "Chebyshev-1"
    # const POINTTYPE = "Chebyshev-2"
    # const POINTTYPE = "Gauss-Legendre"
    # const POINTTYPE = "Gauss-Lobatto"

    const MAXITER = 2000

    #! Data generation
    const SAVEDATA = false # take some time

    #! Ploting constants
    const ADDTITLE = false
    const SAVEPLOT = true
    const FONTSIZE = 20

    "Global parameters"
end

# ╔═╡ 252907f1-b344-4466-8f1a-8ff7ea94d7f7
function green_kernel(k, x, y)
    nxy = norm(x - y)
    return exp(im * k * nxy) / nxy
end

# ╔═╡ f40a28d9-7577-4632-94b7-a617d91c1184
function reference_p(
    fct::Function; sas::src.SelfAffineSet{D,T,N}, nb_pts_cbt::Int
) where {D,T,N}
    cbt1 = src.compute_cubature(sas, "Chebyshev-1", nb_pts_cbt; maxiter=MAXITER)
    cbt2 = src.compute_cubature(sas, "Gauss-Legendre", nb_pts_cbt ÷ 2; maxiter=MAXITER)

    r1, r2 = cbt1(fct), cbt2(fct)

    if isapprox(r1, 0)
        return (r1, abs(r1 - r2))
    end
    return (r1, abs(r2 / r1 - 1))
end

# ╔═╡ c9299045-5429-43ee-a3d1-4a1817f09512
function reference_h(
    fct::Function; sas::src.SelfAffineSet{D,T,N}, nb_pts_cbt::Int, f_diam::T
) where {D,T,N}
    hcbt1 = src.HCubature(
        src.compute_cubature(sas, "Chebyshev-1", nb_pts_cbt; maxiter=MAXITER), sas
    )
    hcbt2 = src.HCubature(
        src.compute_cubature(sas, "Gauss-Legendre", nb_pts_cbt ÷ 2; maxiter=MAXITER), sas
    )

    d = 2 * sas.bounding_ball.radius * f_diam
    while src.diameter(hcbt1) > d
        src.refine!(hcbt1, sas)
        src.refine!(hcbt2, sas)
    end
    r1, r2 = hcbt1(fct), hcbt2(fct)

    if isapprox(r1, 0)
        return (r1, abs(r1 - r2))
    end
    return (r1, abs(r2 / r1 - 1))
end

# ╔═╡ 131d23e1-6451-419f-afce-edd809c7fc67
function sequence_p_version(
    fct::Function; sas::src.SelfAffineSet{D,T,N}, nb_pts_max::Int
) where {D,T,N}
    cbt = src.compute_cubature(sas, POINTTYPE, 2; maxiter=MAXITER)

    nb_pts = [length(cbt)]
    values = [cbt(fct)]

    p = 3
    while length(cbt) < nb_pts_max
        cbt = src.compute_cubature(sas, POINTTYPE, p; maxiter=MAXITER)

        push!(nb_pts, length(cbt))
        push!(values, cbt(fct))

        p += 1
    end

    return (nb_pts, [2 * sas.bounding_ball.radius], values)
end

# ╔═╡ ab675bc7-ecb3-4e39-9350-e0acb5f4be8c
function sequence_h_version(
    fct::Function; sas::src.SelfAffineSet{D,T,N}, cbt::src.Cubature{D,T}, f_diam::T
) where {D,T,N}
    hcbt = src.HCubature(cbt, sas)

    nb_pts = [length(hcbt)]
    mesh_size = [src.diameter(hcbt)]
    values = [hcbt(fct)]

    d = 2 * sas.bounding_ball.radius * f_diam
    while src.diameter(hcbt) > d
        src.refine!(hcbt, sas)

        push!(nb_pts, length(hcbt))
        push!(mesh_size, src.diameter(hcbt))
        push!(values, hcbt(fct))
    end

    return (nb_pts, mesh_size, values)
end

# ╔═╡ 9bbf2b35-a4b1-4943-b519-5149a2be2f1f
function result_to_dict(result)
    return OrderedDict(
        "nb_pts" => result[1],
        "mesh-size" => result[2],
        "values-real" => real.(result[3]),
        "values-imag" => imag.(result[3]),
    )
end

# ╔═╡ b8fdbae2-445d-444f-8561-3fc953a84983
function _save_data(;
    sas::src.SelfAffineSet{D,T,N},
    k::Real,
    x0::SVector{D,T},
    Np::Int,
    f_diam::T,
    suffix::String="",
) where {D,T,N}
    fct = x -> green_kernel(k, x, x0)

    data = OrderedDict("k" => k, "x0" => x0, "point-type" => POINTTYPE)

    result, precision = reference_h(fct; sas=sas, nb_pts_cbt=15, f_diam=f_diam / 2)
    data["reference"] = OrderedDict(
        "value-real" => real(result), "value-imag" => imag(result), "precision" => precision
    )

    data["barycenter-rule"] = result_to_dict(
        sequence_h_version(fct; sas=sas, cbt=src.barycenter_rule(sas), f_diam=f_diam)
    )
    for (i, p) in enumerate(2:6)
        data["h-version-Q$i"] = result_to_dict(
            sequence_h_version(
                fct;
                sas=sas,
                cbt=src.compute_cubature(sas, POINTTYPE, p; maxiter=MAXITER),
                f_diam=f_diam,
            ),
        )
    end
    data["p-version"] = result_to_dict(sequence_p_version(fct; sas=sas, nb_pts_max=Np))

    open("../data-convergences/$(sas.name)$suffix.toml", "w") do io
        TOML.print(io, data)
    end

    return nothing
end

# ╔═╡ ef5ef526-11f9-4a6b-bce1-d21870fea3c1
if SAVEDATA
    for (sas, suffix) in [
        (src.vicsek_2d(1 / 3), ""),
        (src.vicsek_2d(1 / 3, 0.4), "-0.4"),
        (src.vicsek_2d(1 / 3, π / 4), "-pio4"),
        (src.sierpinski_triangle(0.5), ""),
        (src.sierpinski_triangle_fat(2), ""),
        (src.koch_snowflake(), ""),
        (src.cantor_dust_non_sym(), ""),
        (src.barnsley_fern(), ""),
    ]
        _save_data(;
            sas=sas,
            k=5.0,
            x0=SVector{2}([0.1, -2.0]),
            Np=750,
            f_diam=1 / 100,
            suffix=suffix,
        )
    end

    _save_data(;
        sas=src.vicsek_3d(1 / 3, true),
        k=5.0,
        x0=SVector{3}([0.5, -2.0, 1.0]),
        Np=3000,
        f_diam=1 / 100,
    )

    for x in [2.0, 1.5, 1.25, 1.0, 0.0]
        _save_data(;
            sas=src.cantor_dust(1 / 3, [-1.0, 1.0], 2),
            k=5.0,
            x0=SVector{2}([0.1, x]),
            Np=750,
            suffix=@sprintf("-%.2f", x),
            f_diam=1 / 100,
        )
    end
end

# ╔═╡ e7504b5e-edd4-49c2-8227-eaf9019bb8ee
function relative_error(result::Number, reference::Number)
    if isapprox(reference, 0)
        return abs(result)
    end
    return abs(result / reference - 1)
end

# ╔═╡ 5d5bb600-9b89-48b3-90bc-4da6e9f48920
function vicsek_2d_pv()
    fig = Figure(; fontsize=FONTSIZE)

    ax_args = Dict(:xscale => log10, :xlabel => L"M", :yscale => log10)
    if ADDTITLE
        ax_args[:title] = L"$p$-version convergence for the 2d Vicsek"
        ax_args[:ylabel] = L"Relative error$$"
    end
    ax = Axis(fig[1, 1]; ax_args...)

    x_min, x_max = typemax(Int), 0
    for (s, m, leg) in [
        ("", :circle, L"\theta = 0"),
        ("-rot-0.4", :cross, L"\theta = 0.4"),
        ("-rot-pio4", :xcross, L"\theta = \pi / 4"),
    ]
        data = TOML.parsefile("../data-convergences/2d-vicsek$s.toml")

        val_ref = data["reference"]["value-real"] + im .* data["reference"]["value-imag"]

        nb_pts = data["p-version"]["nb_pts"]
        val = data["p-version"]["values-real"] .+ im .* data["p-version"]["values-imag"]

        scatterlines!(
            ax, nb_pts, relative_error.(val, val_ref); marker=m, linestyle=:dash, label=leg
        )

        x_min = min(x_min, nb_pts[1])
        x_max = max(x_max, nb_pts[end])
    end

    limits!(ax, (x_min / 1.1, x_max * 1.1), (1e-16, 10))
    axislegend(ax; position=:lb)

    if SAVEPLOT
        save("2d-vicsek-pv.pdf", fig)
    end

    return fig
end

# ╔═╡ 71c6ed36-30c4-45d9-8cbb-d0d3c7212f62
vicsek_2d_pv()

# ╔═╡ 15ecd841-9e9e-4291-900c-da1f6886a005
function vicsek_2d_hv()
    fig = Figure(; fontsize=FONTSIZE)

    ax_args = Dict(:xscale => log10, :xlabel => L"mesh size$$", :yscale => log10)
    if ADDTITLE
        ax_args[:title] = L"$h$-version convergence for the 2d Vicsek"
        ax_args[:ylabel] = L"Relative error$$"
    end
    ax = Axis(fig[1, 1]; ax_args...)

    h = [2.5e-2, 1.5e-1]
    for (k, y) in [(1, 3e-4), (3, 5e-9), (5, 5e-14)]
        w = y .* (h ./ h[1]) .^ (k + 1)
        lines!(ax, h, w; color=:black)
        p = k + 1
        text!(ax, √prod(h), √prod(w) * 2; text=L"h^%$p")
    end

    for (i, (s, t)) in
        enumerate([("", "0"), ("-rot-0.4", "0.4"), ("-rot-pio4", "\\pi / 4")])
        data = TOML.parsefile("../data-convergences/2d-vicsek$s.toml")

        val_ref = data["reference"]["value-real"] + im .* data["reference"]["value-imag"]

        for (k, ls, mk) in [(1, :dashdot, :circle), (3, :dash, :cross), (5, :dot, :xcross)]
            name = "h-version-Q$k"
            mesh_size = data[name]["mesh-size"]
            val = data[name]["values-real"] .+ im .* data[name]["values-imag"]

            scatterlines!(
                ax,
                mesh_size,
                relative_error.(val, val_ref);
                color=i,
                colormap=:tab10,
                colorrange=(1, 10),
                marker=mk,
                linestyle=ls,
                label=L"$\mathbb{Q}_%$k$, $\theta = %$t $",
            )
        end
    end

    limits!(ax, (1e-2, 3.5), (1e-16, 10))
    axislegend(ax; position=:rb, backgroundcolor=(:white, 0))

    if SAVEPLOT
        save("2d-vicsek-hv.pdf", fig)
    end

    return fig
end

# ╔═╡ 9356d11c-bf4c-4960-90e6-b787ed55a3e4
vicsek_2d_hv()

# ╔═╡ 53324fc5-0675-48be-a9f0-1da0e8ef5e10
function other_example_pv()
    fig = Figure(; fontsize=FONTSIZE)

    ax_args = Dict(:xlabel => L"M", :yscale => log10)
    if ADDTITLE
        ax_args[:ylabel] = L"Relative error$$"
    end
    ax = Axis(fig[1, 1]; ax_args...)

    for (i, (name, mk, leg)) in enumerate([
        ("2d-sierpinski-triangle-fat", :circle, "fat Sierpinski tri."),
        ("2d-koch-snowflake", :cross, "Koch snowflake"),
        ("2d-cantor-non-sym", :xcross, "Cantor dust"),
    ])
        data = TOML.parsefile("../data-convergences/$name.toml")

        val_ref = data["reference"]["value-real"] + im .* data["reference"]["value-imag"]

        nb_pts = data["p-version"]["nb_pts"]
        val = data["p-version"]["values-real"] .+ im .* data["p-version"]["values-imag"]

        scatterlines!(
            ax,
            nb_pts,
            relative_error.(val, val_ref);
            color=i,
            colormap=:tab10,
            colorrange=(1, 10),
            marker=mk,
            linestyle=:dash,
            label=L"%$leg$$",
        )
    end

    limits!(ax, (0, 800), (1e-16, 10))
    axislegend(ax; position=:rt, backgroundcolor=(:white, 0))

    if SAVEPLOT
        save("2d-other-example-pv.pdf", fig)
    end

    return fig
end

# ╔═╡ 741e4d28-56dd-4035-b543-3114173be543
other_example_pv()

# ╔═╡ 8e335283-2a0b-4677-9632-7bc17e583f4f
function other_example_hv()
    fig = Figure(; fontsize=FONTSIZE)

    ax_args = Dict(
        :xscale => log10,
        :xlabel => L"mesh size$$",
        :yscale => log10,
        :ylabel => L"Relative error$$",
    )
    if ADDTITLE
        ax_args[:title] = L"$h$-version convergence for %$name"
    end
    ax = Axis(fig[1, 1]; ax_args...)

    h = [2.5e-2, 1.5e-1]
    for (k, y) in [(1, 2e-3), (3, 3e-8), (5, 5e-13)]
        w = y .* (h ./ h[1]) .^ (k + 1)
        lines!(ax, h, w; color=:black)
        p = k + 1
        text!(ax, √prod(h), √prod(w) * 2; text=L"h^%$p")
    end

    for (i, (filename, leg)) in enumerate([
        ("2d-sierpinski-triangle-fat", "fat Sierpinski tri."),
        ("2d-koch-snowflake", "Koch snowflake"),
        ("2d-cantor-non-sym", "Cantor dust"),
    ])
        data = TOML.parsefile("../data-convergences/$filename.toml")

        val_ref = data["reference"]["value-real"] + im .* data["reference"]["value-imag"]

        for (k, mk) in [(1, :circle), (3, :cross), (5, :xcross)]
            name = "h-version-Q$k"
            mesh_size = data[name]["mesh-size"]
            val = data[name]["values-real"] .+ im .* data[name]["values-imag"]

            scatterlines!(
                ax,
                mesh_size,
                relative_error.(val, val_ref);
                color=i,
                colormap=:tab10,
                colorrange=(1, 10),
                marker=mk,
                linestyle=:dash,
                label=L"$\mathbb{Q}_%$k$, %$leg",
            )
        end
    end

    limits!(ax, (1e-2, 10), (1e-16, 10))
    axislegend(ax; position=:rb, backgroundcolor=(:white, 0))

    if SAVEPLOT
        save("2d-other-example-hv.pdf", fig)
    end

    return fig
end

# ╔═╡ 43eb7cbf-c71e-4b7f-aa91-3fd5348de479
other_example_hv()

# ╔═╡ 44c81e00-8cd3-409c-9d0f-280482d792f4
function vicsek_3d_pv()
    fig = Figure(; fontsize=FONTSIZE)

    ax_args = Dict(:xlabel => L"M", :yscale => log10)
    if ADDTITLE
        ax_args[:title] = L"$p$-version convergence for the 3d Vicsek"
        ax_args[:ylabel] = L"Relative error$$"
    end
    ax = Axis(fig[1, 1]; ax_args...)

    data = TOML.parsefile("../data-convergences/3d-vicsek-rot.toml")

    val_ref = data["reference"]["value-real"] + im .* data["reference"]["value-imag"]

    nb_pts = data["p-version"]["nb_pts"]
    val = data["p-version"]["values-real"] .+ im .* data["p-version"]["values-imag"]

    scatterlines!(
        ax, nb_pts, relative_error.(val, val_ref); marker=:xcross, linestyle=:dash
    )

    limits!(ax, (0, 3500), (1e-16, 10))

    if SAVEPLOT
        save("3d-vicsek-rot-pv.pdf", fig)
    end

    return fig
end

# ╔═╡ e5f91725-d93c-4a0d-b60a-84c1aa98c970
vicsek_3d_pv()

# ╔═╡ 1c6e496b-c5e0-46c5-bd6e-7103fabba89a
function vicsek_3d_hv()
    fig = Figure(; fontsize=FONTSIZE)

    ax_args = Dict(:xscale => log10, :xlabel => L"mesh size$$", :yscale => log10)
    if ADDTITLE
        ax_args[:title] = L"$h$-version convergence for %$name"
        ax_args[:ylabel] = L"Relative error$$"
    end
    ax = Axis(fig[1, 1]; ax_args...)

    h = [2.5e-2, 1.5e-1]
    for (k, y) in [(1, 2e-4), (3, 1e-8), (5, 3e-13)]
        w = y .* (h ./ h[1]) .^ (k + 1)
        lines!(ax, h, w; color=:black)
        p = k + 1
        text!(ax, √prod(h), √prod(w) * 2; text=L"h^%$p")
    end

    data = TOML.parsefile("../data-convergences/3d-vicsek-rot.toml")

    val_ref = data["reference"]["value-real"] + im .* data["reference"]["value-imag"]

    for (i, k) in enumerate([1, 3, 5])
        name = "h-version-Q$k"
        mesh_size = data[name]["mesh-size"]
        val = data[name]["values-real"] .+ im .* data[name]["values-imag"]

        scatterlines!(
            ax,
            mesh_size,
            relative_error.(val, val_ref);
            color=i,
            colormap=:tab10,
            colorrange=(1, 10),
            marker=:xcross,
            linestyle=:dash,
            label=L"$\mathbb{Q}_%$k$",
        )
    end

    limits!(ax, (1e-2, 4), (1e-16, 10))
    axislegend(ax; position=:rb)

    if SAVEPLOT
        save("3d-vicsek-rot-hv.pdf", fig)
    end

    return fig
end

# ╔═╡ 6dbb183a-2f10-4a58-89de-cf4d0f98eec6
vicsek_3d_hv()

# ╔═╡ 761fb9fd-d4df-40fc-a30f-50a40695a129
function barnsley_fern_pv()
    fig = Figure(; fontsize=FONTSIZE)

    ax_args = Dict(:xlabel => L"M", :yscale => log10)
    if ADDTITLE
        ax_args[:title] = L"$p$-version convergence for the 3d Vicsek"
        ax_args[:ylabel] = L"Relative error$$"
    end
    ax = Axis(fig[1, 1]; ax_args...)

    data = TOML.parsefile("../data-convergences/2d-barnsley-fern.toml")

    val_ref = data["reference"]["value-real"] + im .* data["reference"]["value-imag"]

    nb_pts = data["p-version"]["nb_pts"]
    val = data["p-version"]["values-real"] .+ im .* data["p-version"]["values-imag"]

    scatterlines!(
        ax, nb_pts, relative_error.(val, val_ref); marker=:xcross, linestyle=:dash
    )

    limits!(ax, (0, 800), (1e-10, 20))

    if SAVEPLOT
        save("2d-barnsley-fern-pv.pdf", fig)
    end

    return fig
end

# ╔═╡ 702de76b-2122-43e0-99c3-7a26a24e77f7
barnsley_fern_pv()

# ╔═╡ 57536c2b-09ef-4c36-abef-3802628a798f
function barnsley_fern_hv()
    fig = Figure(; fontsize=FONTSIZE)

    ax_args = Dict(:xscale => log10, :xlabel => L"mesh size$$", :yscale => log10)
    if ADDTITLE
        ax_args[:title] = L"$h$-version convergence for %$name"
        ax_args[:ylabel] = L"Relative error$$"
    end
    ax = Axis(fig[1, 1]; ax_args...)

    h = [1e-1, 4e-1]
    for (k, y) in [(1, 1e-2), (3, 5e-6), (5, 2e-9)]
        w = y .* (h ./ h[1]) .^ (k + 1)
        lines!(ax, h, w; color=:black)
        p = k + 1
        text!(ax, √prod(h), √prod(w) * 2; text=L"h^%$p")
    end

    data = TOML.parsefile("../data-convergences/2d-barnsley-fern.toml")

    val_ref = data["reference"]["value-real"] + im .* data["reference"]["value-imag"]

    for (i, k) in enumerate([1, 3, 5])
        name = "h-version-Q$k"
        mesh_size = data[name]["mesh-size"]
        val = data[name]["values-real"] .+ im .* data[name]["values-imag"]

        scatterlines!(
            ax,
            mesh_size,
            relative_error.(val, val_ref);
            color=i,
            colormap=:tab10,
            colorrange=(1, 10),
            marker=:xcross,
            linestyle=:dash,
            label=L"$\mathbb{Q}_%$k$",
        )
    end

    limits!(ax, (8e-2, 15), (1e-10, 20))
    axislegend(ax; position=:rb)

    if SAVEPLOT
        save("2d-barnsley-fern-hv.pdf", fig)
    end

    return fig
end

# ╔═╡ 8decb395-ecdb-4665-b2b3-9cb16f8c6c75
barnsley_fern_hv()

# ╔═╡ 2f67fe32-1768-4013-bc81-52350290242a
function cantor_dust_sing_pv()
    fig = Figure(; fontsize=FONTSIZE)

    ax_args = Dict(:xlabel => L"number of cubature points$$", :yscale => log10)
    if ADDTITLE
        ax_args[:title] = L"$p$-version convergence for nonsymmetric Cantor dust"
        ax_args[:ylabel] = L"Relative error$$"
    end
    ax = Axis(fig[1, 1]; ax_args...)

    for (x, c, mk) in [
        (-2.0, 1, :circle),
        (-1.5, 2, :cross),
        (-1.25, 3, :xcross),
        (-1.0, 4, :utriangle),
        (0.0, 5, :diamond),
    ]
        filename = "2d-cantor-dust-$(@sprintf("%.2f", abs(x)))"
        data = TOML.parsefile("../data-convergences/$filename.toml")

        val_ref = data["reference"]["value-real"] + im .* data["reference"]["value-imag"]

        nb_pts = data["p-version"]["nb_pts"]
        val = data["p-version"]["values-real"] .+ im .* data["p-version"]["values-imag"]

        scatterlines!(
            ax,
            nb_pts,
            relative_error.(val, val_ref);
            color=c,
            colormap=:tab10,
            colorrange=(1, 10),
            marker=mk,
            linestyle=:dash,
            label=L"y = %$x",
        )
    end

    limits!(ax, (0, 800), (1e-16, 10))
    axislegend(ax; position=:lb)

    if SAVEPLOT
        save("2d-cantor-dust-sing-pv.pdf", fig)
    end

    return fig
end

# ╔═╡ 13371d5a-fb51-4010-b139-49617b7d314d
cantor_dust_sing_pv()

# ╔═╡ f19cafcc-b3e2-4998-b08b-cc57606f8fba
function cantor_dust_sing_hv()
    fig = Figure(; fontsize=FONTSIZE)

    ax_args = Dict(:xscale => log10, :xlabel => L"mesh size$$", :yscale => log10)
    if ADDTITLE
        ax_args[:title] = L"$p$-version convergence for nonsymmetric Cantor dust"
        ax_args[:ylabel] = L"Relative error$$"
    end
    ax = Axis(fig[1, 1]; ax_args...)

    h = [2.5e-2, 1.5e-1]
    for (k, y) in [(1, 4e-4), (3, 8e-8), (5, 3e-12)]
        w = y .* (h ./ h[1]) .^ (k + 1)
        lines!(ax, h, w; color=:black)
        p = k + 1
        text!(ax, √prod(h), √prod(w) * 2; text=L"h^%$p")
    end

    for (i, x) in enumerate([-2.0, -1.5, -1.25, -1.0, 0.0])
        filename = "2d-cantor-dust-$(@sprintf("%.2f", abs(x)))"
        data = TOML.parsefile("../data-convergences/$filename.toml")

        val_ref = data["reference"]["value-real"] + im .* data["reference"]["value-imag"]

        for (k, mk, ls) in [(1, :circle, :dashdot), (3, :cross, :dash), (5, :xcross, :dot)]
            name = "h-version-Q$k"
            mesh_size = data[name]["mesh-size"]
            val = data[name]["values-real"] .+ im .* data[name]["values-imag"]

            args = Dict(
                :color => k,
                :colormap => :tab10,
                :colorrange => (1, 10),
                :marker => mk,
                :linestyle => ls,
            )
            if i == 1
                args[:label] = L"$\mathbb{Q}_%$k$"
            end
            scatterlines!(ax, mesh_size, relative_error.(val, val_ref); args...)
        end
    end

    limits!(ax, (1e-2, 3.5), (1e-16, 10))
    axislegend(ax; position=:rb)

    if SAVEPLOT
        save("2d-cantor-dust-sing-hv.pdf", fig)
    end

    return fig
end

# ╔═╡ 2b9c98a7-218a-4770-ad66-1c7f0cc450c5
cantor_dust_sing_hv()

# ╔═╡ a0b5d444-1dec-4d63-ab09-bd3ee8e02bc3
function _plot_pv(name::String)
    fig = Figure(; fontsize=FONTSIZE)

    ax_args = Dict(
        :xlabel => L"number of cubature points$$",
        :yscale => log10,
        :ylabel => L"Relative error$$",
    )
    if ADDTITLE
        ax_args[:title] = L"$p$-version convergence for %$name"
    end
    ax = Axis(fig[1, 1]; ax_args...)

    data = TOML.parsefile("../data-convergences/$name.toml")

    val_ref = data["reference"]["value-real"] + im .* data["reference"]["value-imag"]

    nb_pts = data["p-version"]["nb_pts"]
    val = data["p-version"]["values-real"] .+ im .* data["p-version"]["values-imag"]

    scatterlines!(
        ax, nb_pts, relative_error.(val, val_ref); marker=:xcross, linestyle=:dash
    )

    return fig
end

# ╔═╡ f281cd8d-b427-4e26-ad91-8538c6e5f744
_plot_pv("2d-sierpinski-triangle")

# ╔═╡ 84dd7bad-7797-45e4-b2e0-f8defab92759
_plot_pv("2d-sierpinski-triangle-fat")

# ╔═╡ 90b65bf1-a164-401c-92b6-5563eb768474
_plot_pv("3d-vicsek-rot")

# ╔═╡ a3b65e8e-14d6-4c0c-a4ef-a29c9695c2bd
_plot_pv("2d-koch-snowflake")

# ╔═╡ 7038e5de-a1cb-42d6-b21b-2d4f77656f62
_plot_pv("2d-cantor-non-sym")

# ╔═╡ a3b8cf1e-eb3f-403e-91ac-8d09f15b8f75
_plot_pv("2d-barnsley-fern")

# ╔═╡ 12535b8a-b083-4704-b889-db0a460120ff
function _plot_hv(name::String)
    fig = Figure(; fontsize=FONTSIZE)

    ax_args = Dict(
        :xscale => log10,
        :xlabel => L"mesh size$$",
        :yscale => log10,
        :ylabel => L"Relative error$$",
    )
    if ADDTITLE
        ax_args[:title] = L"$h$-version convergence for %$name"
    end
    ax = Axis(fig[1, 1]; ax_args...)

    data = TOML.parsefile("../data-convergences/$name.toml")

    val_ref = data["reference"]["value-real"] + im .* data["reference"]["value-imag"]

    for k in 1:5
        name = "h-version-Q$k"
        mesh_size = data[name]["mesh-size"]
        val = data[name]["values-real"] .+ im .* data[name]["values-imag"]

        scatterlines!(
            ax,
            mesh_size,
            relative_error.(val, val_ref);
            color=k,
            colormap=:tab10,
            colorrange=(1, 10),
            marker=:xcross,
            linestyle=:dash,
            label=L"$Q%$k$",
        )
    end

    axislegend(ax; position=:rb)

    return fig
end

# ╔═╡ 50c9daf2-7ebc-434f-a582-06c1cf6d0e4f
_plot_hv("2d-sierpinski-triangle")

# ╔═╡ 95cbcb72-0fda-48b2-8355-081f4c38f5b6
_plot_hv("2d-sierpinski-triangle-fat")

# ╔═╡ 2876ac50-8d20-4e82-8684-f9f626952fa0
_plot_hv("3d-vicsek-rot")

# ╔═╡ 732a6a47-ede0-432f-bbf6-dbf8cd11f40d
_plot_hv("2d-koch-snowflake")

# ╔═╡ e4d4001a-dc17-413e-b118-e0317c170a0d
_plot_hv("2d-cantor-non-sym")

# ╔═╡ c3af3f2b-b49a-45eb-9db3-37d674f1a48c
_plot_hv("2d-barnsley-fern")

# ╔═╡ Cell order:
# ╠═fe657a18-63c0-11ef-0959-15170cb4d30a
# ╠═b5d45ad9-3996-401f-96e5-78aba17e8686
# ╠═ef5ef526-11f9-4a6b-bce1-d21870fea3c1
# ╠═5d5bb600-9b89-48b3-90bc-4da6e9f48920
# ╠═71c6ed36-30c4-45d9-8cbb-d0d3c7212f62
# ╠═15ecd841-9e9e-4291-900c-da1f6886a005
# ╠═9356d11c-bf4c-4960-90e6-b787ed55a3e4
# ╠═53324fc5-0675-48be-a9f0-1da0e8ef5e10
# ╠═741e4d28-56dd-4035-b543-3114173be543
# ╠═8e335283-2a0b-4677-9632-7bc17e583f4f
# ╠═43eb7cbf-c71e-4b7f-aa91-3fd5348de479
# ╠═44c81e00-8cd3-409c-9d0f-280482d792f4
# ╠═e5f91725-d93c-4a0d-b60a-84c1aa98c970
# ╠═1c6e496b-c5e0-46c5-bd6e-7103fabba89a
# ╠═6dbb183a-2f10-4a58-89de-cf4d0f98eec6
# ╠═761fb9fd-d4df-40fc-a30f-50a40695a129
# ╠═702de76b-2122-43e0-99c3-7a26a24e77f7
# ╠═57536c2b-09ef-4c36-abef-3802628a798f
# ╠═8decb395-ecdb-4665-b2b3-9cb16f8c6c75
# ╠═2f67fe32-1768-4013-bc81-52350290242a
# ╠═13371d5a-fb51-4010-b139-49617b7d314d
# ╠═f19cafcc-b3e2-4998-b08b-cc57606f8fba
# ╠═2b9c98a7-218a-4770-ad66-1c7f0cc450c5
# ╠═a0b5d444-1dec-4d63-ab09-bd3ee8e02bc3
# ╠═12535b8a-b083-4704-b889-db0a460120ff
# ╠═f281cd8d-b427-4e26-ad91-8538c6e5f744
# ╠═50c9daf2-7ebc-434f-a582-06c1cf6d0e4f
# ╠═84dd7bad-7797-45e4-b2e0-f8defab92759
# ╠═95cbcb72-0fda-48b2-8355-081f4c38f5b6
# ╠═90b65bf1-a164-401c-92b6-5563eb768474
# ╠═2876ac50-8d20-4e82-8684-f9f626952fa0
# ╠═a3b65e8e-14d6-4c0c-a4ef-a29c9695c2bd
# ╠═732a6a47-ede0-432f-bbf6-dbf8cd11f40d
# ╠═7038e5de-a1cb-42d6-b21b-2d4f77656f62
# ╠═e4d4001a-dc17-413e-b118-e0317c170a0d
# ╠═a3b8cf1e-eb3f-403e-91ac-8d09f15b8f75
# ╠═c3af3f2b-b49a-45eb-9db3-37d674f1a48c
# ╠═b8fdbae2-445d-444f-8561-3fc953a84983
# ╠═252907f1-b344-4466-8f1a-8ff7ea94d7f7
# ╠═f40a28d9-7577-4632-94b7-a617d91c1184
# ╠═c9299045-5429-43ee-a3d1-4a1817f09512
# ╠═131d23e1-6451-419f-afce-edd809c7fc67
# ╠═ab675bc7-ecb3-4e39-9350-e0acb5f4be8c
# ╠═9bbf2b35-a4b1-4943-b519-5149a2be2f1f
# ╠═e7504b5e-edd4-49c2-8227-eaf9019bb8ee
