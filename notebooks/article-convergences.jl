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
    const ADDTITLE = true
    const SAVEPLOT = false
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
    hcbt1 = src.HCubature(src.compute_cubature(sas, "Chebyshev-1", nb_pts_cbt), sas)
    hcbt2 = src.HCubature(src.compute_cubature(sas, "Gauss-Legendre", nb_pts_cbt ÷ 2), sas)

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
                fct; sas=sas, cbt=src.compute_cubature(sas, POINTTYPE, p), f_diam=f_diam
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
    for sas in [
        src.vicsek_2d(1 / 3),
        src.vicsek_2d(1 / 3, π / 4),
        src.sierpinski_triangle(0.5),
        src.sierpinski_triangle_fat(2),
        src.koch_snowflake(),
        src.cantor_dust_non_sym(),
        src.barnsley_fern(),
    ]
        _save_data(; sas=sas, k=5.0, x0=SVector{2}([0.5, -2.0]), Np=750, f_diam=1 / 100)
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
            x0=SVector{2}([x, 0.1]),
            Np=626,
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
function vicsek_pv()
    fig = Figure(; fontsize=FONTSIZE)

    ax_args = Dict(
        :xscale => log10, :xlabel => L"M", :yscale => log10, :ylabel => L"Relative error$$"
    )
    if ADDTITLE
        ax_args[:title] = L"$p$-version convergence for the Vicsek"
    end
    ax = Axis(fig[1, 1]; ax_args...)

    x_min, x_max = typemax(Int), 0
    for (s, m, l) in [("", :cross, "without"), ("-rot", :xcross, "with")]
        data = TOML.parsefile("../data-convergences/2d-vicsek$s.toml")

        val_ref = data["reference"]["value-real"] + im .* data["reference"]["value-imag"]

        nb_pts = data["p-version"]["nb_pts"]
        val = data["p-version"]["values-real"] .+ im .* data["p-version"]["values-imag"]

        scatterlines!(
            ax,
            nb_pts,
            relative_error.(val, val_ref);
            marker=m,
            linestyle=:dash,
            label=L"%$l rotation$$",
        )

        x_min = min(x_min, nb_pts[1])
        x_max = max(x_max, nb_pts[end])
    end

    limits!(ax, (x_min / 1.1, x_max * 1.1), (1e-16, 10))
    axislegend(ax; position=:rt)

    if SAVEPLOT
        save("2d-vicsek-pv.pdf", fig)
    end

    return fig
end

# ╔═╡ 71c6ed36-30c4-45d9-8cbb-d0d3c7212f62
vicsek_pv()

# ╔═╡ 15ecd841-9e9e-4291-900c-da1f6886a005
function vicsek_hv()
    fig = Figure(; fontsize=FONTSIZE)

    ax_args = Dict(
        :xscale => log10, :xlabel => L"M", :yscale => log10, :ylabel => L"Relative error$$"
    )
    if ADDTITLE
        ax_args[:title] = L"$h$-version convergence for the Vicsek"
    end
    ax = Axis(fig[1, 1]; ax_args...)

    x_min, x_max = typemax(Int), 0
    for (s, mk, ls, leg) in [("", :cross, :dash, ""), ("-rot", :xcross, :dot, "rot")]
        data = TOML.parsefile("../data-convergences/2d-vicsek$s.toml")

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

            x_min = min(x_min, mesh_size[end])
            x_max = max(x_max, mesh_size[1])
        end
    end

    limits!(ax, (x_min / 1.1, x_max * 1.1), (1e-16, 10))
    axislegend(ax; position=:rb)

    if SAVEPLOT
        save("2d-vicsek-hv.pdf", fig)
    end

    return fig
end

# ╔═╡ 9356d11c-bf4c-4960-90e6-b787ed55a3e4
vicsek_hv()

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

    x_min, x_max = typemax(Int), 0
    data = TOML.parsefile("../data-convergences/$name.toml")

    val_ref = data["reference"]["value-real"] + im .* data["reference"]["value-imag"]

    nb_pts = data["p-version"]["nb_pts"]
    val = data["p-version"]["values-real"] .+ im .* data["p-version"]["values-imag"]

    scatterlines!(
        ax, nb_pts, relative_error.(val, val_ref); marker=:xcross, linestyle=:dash
    )

    x_min = min(x_min, nb_pts[1])
    x_max = max(x_max, nb_pts[end])

    # limits!(ax, (x_min / 1.1, x_max * 1.1), (1e-16, 10))

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

    x_min, x_max = typemax(Int), 0
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

        x_min = min(x_min, mesh_size[end])
        x_max = max(x_max, mesh_size[1])
    end

    # limits!(ax, (x_min / 1.1, x_max * 1.1), (1e-16, 10))
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

# ╔═╡ bf9d9043-4c1c-4b82-926d-a3ee7aa0574f
function cantor_dust_sing_pv()
    fig = Figure(; fontsize=FONTSIZE)

    ax_args = Dict(
        :xlabel => L"number of cubature points$$",
        :yscale => log10,
        :ylabel => L"Relative error$$",
    )
    if ADDTITLE
        ax_args[:title] = L"$p$-version convergence for nonsymmetric Cantor dust"
    end
    ax = Axis(fig[1, 1]; ax_args...)

    x_min, x_max = typemax(Int), 0
    for (i, x) in enumerate([2.0, 1.5, 1.25, 1.0, 0.0])
        filename = "2d-cantor-dust-$(@sprintf("%.2f", x))"
        data = TOML.parsefile("../data-convergences/$filename.toml")

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
            marker=:xcross,
            linestyle=:dash,
            label=L"x = %$x",
        )

        x_min = min(x_min, nb_pts[1])
        x_max = max(x_max, nb_pts[end])
    end

    limits!(ax, (x_min / 1.1, x_max * 1.1), (1e-16, 10))
    axislegend(ax; position=:lb)

    return fig
end

# ╔═╡ 77fe701f-f2ce-408b-a785-bb9acee6a54a
cantor_dust_sing_pv()

# ╔═╡ d1814922-76e9-4469-93f9-42c2a27bf98b
function cantor_dust_sing_hv(k::Int)
    fig = Figure(; fontsize=FONTSIZE)

    ax_args = Dict(
        :xscale => log10, :xlabel => L"M", :yscale => log10, :ylabel => L"Relative error$$"
    )
    if ADDTITLE
        ax_args[:title] = L"$p$-version convergence for nonsymmetric Cantor dust"
    end
    ax = Axis(fig[1, 1]; ax_args...)

    x_min, x_max = typemax(Int), 0
    for (i, x) in enumerate([2.0, 1.5, 1.25, 1.0, 0.0])
        filename = "2d-cantor-dust-$(@sprintf("%.2f", x))"
        data = TOML.parsefile("../data-convergences/$filename.toml")

        val_ref = data["reference"]["value-real"] + im .* data["reference"]["value-imag"]

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

        x_min = min(x_min, mesh_size[end])
        x_max = max(x_max, mesh_size[1])
    end

    limits!(ax, (x_min / 1.1, x_max * 1.1), (1e-16, 10))
    axislegend(ax; position=:rb)

    return fig
end

# ╔═╡ fef04d01-a9bf-4cef-9bc3-a3cacd92d49d
cantor_dust_sing_hv(4)

# ╔═╡ Cell order:
# ╠═fe657a18-63c0-11ef-0959-15170cb4d30a
# ╠═b5d45ad9-3996-401f-96e5-78aba17e8686
# ╠═ef5ef526-11f9-4a6b-bce1-d21870fea3c1
# ╠═5d5bb600-9b89-48b3-90bc-4da6e9f48920
# ╠═71c6ed36-30c4-45d9-8cbb-d0d3c7212f62
# ╠═15ecd841-9e9e-4291-900c-da1f6886a005
# ╠═9356d11c-bf4c-4960-90e6-b787ed55a3e4
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
# ╠═bf9d9043-4c1c-4b82-926d-a3ee7aa0574f
# ╠═77fe701f-f2ce-408b-a785-bb9acee6a54a
# ╠═d1814922-76e9-4469-93f9-42c2a27bf98b
# ╠═fef04d01-a9bf-4cef-9bc3-a3cacd92d49d
# ╠═b8fdbae2-445d-444f-8561-3fc953a84983
# ╠═252907f1-b344-4466-8f1a-8ff7ea94d7f7
# ╠═f40a28d9-7577-4632-94b7-a617d91c1184
# ╠═c9299045-5429-43ee-a3d1-4a1817f09512
# ╠═131d23e1-6451-419f-afce-edd809c7fc67
# ╠═ab675bc7-ecb3-4e39-9350-e0acb5f4be8c
# ╠═9bbf2b35-a4b1-4943-b519-5149a2be2f1f
# ╠═e7504b5e-edd4-49c2-8227-eaf9019bb8ee
