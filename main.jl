using Printf, LinearAlgebra, StaticArrays, DataStructures, TOML, CairoMakie

import IFSCubature as src

const POINTTYPE = "Chebyshev-1"

const MAXITER = 2048

#! Ploting constants
const ADDTITLE = true
const SAVEPLOT = false
const FONTSIZE = 20

function green_kernel(k, x, y)
    nxy = norm(x - y)
    return exp(im * k * nxy) / nxy
end

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

function result_to_dict(result)
    return OrderedDict(
        "nb_pts" => result[1],
        "mesh-size" => result[2],
        "values-real" => real.(result[3]),
        "values-imag" => imag.(result[3]),
    )
end

function _save_data(;
    sas::src.SelfAffineSet{D,T,N}, k::Real, x0::SVector{D,T}, Np::Int, suffix::String=""
) where {D,T,N}
    fct = x -> green_kernel(k, x, x0)

    data = OrderedDict("k" => k, "x0" => x0, "point-type" => POINTTYPE)

    result, precision = reference_h(fct; sas=sas, nb_pts_cbt=15)
    data["reference"] = OrderedDict(
        "value-real" => real(result), "value-imag" => imag(result), "precision" => precision
    )

    data["barycenter-rule"] = result_to_dict(
        sequence_h_version(fct; sas=sas, cbt=src.barycenter_rule(sas))
    )
    for (i, p) in enumerate(2:6)
        data["h-version-Q$i"] = result_to_dict(
            sequence_h_version(fct; sas=sas, cbt=src.compute_cubature(sas, POINTTYPE, p))
        )
    end
    data["p-version"] = result_to_dict(sequence_p_version(fct; sas=sas, nb_pts_max=Np))

    open("./data-convergences/$(sas.name)$suffix.toml", "w") do io
        TOML.print(io, data)
    end

    return nothing
end

function relative_error(result::Number, reference::Number)
    if isapprox(reference, 0)
        return abs(result)
    end
    return abs(result / reference - 1)
end

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
    data = TOML.parsefile("./data-convergences/$name.toml")

    val_ref = data["reference"]["value-real"] + im .* data["reference"]["value-imag"]

    nb_pts = data["p-version"]["nb_pts"]
    val = data["p-version"]["values-real"] .+ im .* data["p-version"]["values-imag"]

    scatterlines!(
        ax, nb_pts, relative_error.(val, val_ref); marker=:xcross, linestyle=:dash
    )

    x_min = min(x_min, nb_pts[1])
    x_max = max(x_max, nb_pts[end])

    limits!(ax, (x_min / 1.1, x_max * 1.1), (1e-16, 10))

    return fig
end

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
    data = TOML.parsefile("./data-convergences/$name.toml")

    val_ref = data["reference"]["value-real"] + im .* data["reference"]["value-imag"]

    for k in 1:5
        name = "h-version-Q$k"
        # nb_pts = data[name]["nb_pts"]
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

function cantor_dust_test()
    ifs = [
        src.contractive_similarity(ρ, c) for (ρ, c) in
        [(0.4, [-1.0, -1.0]), (0.3, [1.0, -1.0]), (0.2, [1.0, 1.0]), (0.1, [-1.0, 1.0])]
    ]

    ball = src.hyper_ball(fill(0.0, 2), √2)
    box = src.hyper_box(fill(0.0, 2), Matrix(Diagonal(fill(1.0, 2))))

    return src.SelfAffineSet(ifs, fill(1 / 4, 4), ball, box, "2d-cantor-dust-test")
end

sas = src.square(-1.0, 1.0)

fig = Figure()
ax = Axis(fig[1, 1]; yscale=log10)

lines!(range(1, 800, 64), exp.(-range(1, 800, 64) .^ 0.5); color=:black)

for y0 in [2.0, 1.5, 1.25, 1.125]
    k = 5.0
    fct = x -> green_kernel(k, x, SVector{2}([0.1, -y0]))

    # _save_data(; sas=sas, k=5.0, x0=SVector{2}([2.0, 0.5]), Np=626)

    # display(_plot_pv(sas.name))
    # display(_plot_hv(sas.name))

    # println("------------------------------------------")
    # v_ref, prec = reference_p(fct; sas=sas, nb_pts_cbt=32)
    # println((v_ref, prec))
    v_ref, prec = reference_h(fct; sas=sas, nb_pts_cbt=15, f_diam=1 / 200)
    println((v_ref, prec))

    # cbt = src.compute_cubature(sas, "Chebyshev-1", 5; maxiter=MAXITER)
    # n, h, v = sequence_h_version(fct; sas=sas, cbt=cbt, f_diam=1 / 100)

    n, h, v = sequence_p_version(fct; sas=sas, nb_pts_max=750)

    scatterlines!(ax, n, abs.(v ./ v_ref .- 1); linestyle=:dash, label="$y0")
end

axislegend(ax; position=:lb)

display(fig)
