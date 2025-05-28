### A Pluto.jl notebook ###
# v0.20.9

using Markdown
using InteractiveUtils

# ╔═╡ a6000660-39a4-11f0-139e-59028a2cdaa8
begin
    using Pkg: Pkg
    Pkg.activate(Base.current_project())
    Pkg.instantiate()

    using Printf, StaticArrays, LinearAlgebra, CairoMakie

    import IFSCubature as src
end

# ╔═╡ 190c44f9-f13c-4d6d-be84-b250e0f3ff61
begin
    #! Type of points
    # const POINTTYPE = "Equispaced-1"
    # const POINTTYPE = "Equispaced-2"
    const POINTTYPE = "Chebyshev-1"
    # const POINTTYPE = "Chebyshev-2"
    # const POINTTYPE = "Gauss-Legendre"
    # const POINTTYPE = "Gauss-Lobatto"

    const MAXITER = 2000

    #! Ploting constants
    const FONTSIZE = 20
    const SAVEFIG = true

    "Global parameters"
end

# ╔═╡ 0a613021-6729-484d-b919-7690a45c6493
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

# ╔═╡ 0a1f6f87-f80c-4aaf-9162-3495b00152c1
function plot_vicsek(xs_list::Vector{Vector{Float64}}, pts_chaos_nb::Int, α_attractor::Real)
    fig = Figure(; size=(600, 400), fontsize=FONTSIZE)
    ax = Axis(fig[1, 1]; aspect=DataAspect(), xlabel=L"x", ylabel=L"y")
    limits!(ax, -2.05, 1.05, -1.05, 1.05)

    poly!(
        ax,
        [SVector(-1, -1), SVector(1, -1), SVector(1, 1), SVector(-1, 1)];
        color=(:black, 0.1),
    )

    sas = src.vicsek_2d(1 / 3)
    plot_chaos_game!(ax, sas, pts_chaos_nb, α_attractor)

    colors, i = Makie.to_colormap(:tab10), 0
    for xs in xs_list
        for x in xs
            i += 1
            scatter!(ax, [x], [-0.1]; color=colors[i], markersize=16)
        end

        if SAVEFIG
            save("2d-vicsek-singular-$i.pdf", fig)
        end
    end

    return fig
end

# ╔═╡ 26308e9f-52b5-465b-b121-96ee8c50aab3
plot_vicsek([[-2.0], [-1.4, -1.2], [-1.0, -0.67]], 200_000, 0.75)

# ╔═╡ 033139d9-6c56-4e10-924f-306c0b2fa627
function green_kernel(k, x, y)
    nxy = norm(x - y)
    return cis(k * nxy) / nxy
end

# ╔═╡ 2c55c152-b2e3-432a-a781-b6045f12e5e2
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

# ╔═╡ 1641cc06-e7f2-4d39-89fe-0446c6edbfe9
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

# ╔═╡ 358da4c0-022b-4232-ad23-9e1e34729463
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

# ╔═╡ 21180f74-48bf-4a85-a818-3127799681d4
function relative_error(result::Number, reference::Number)
    if isapprox(reference, 0)
        return abs(result)
    end
    return abs(result / reference - 1)
end

# ╔═╡ fe3d2924-9fd2-4fd9-b0e3-6522e1d95edf
function plot_vicsek_pv_sing(xs_list::Vector{Vector{Float64}})
    fig = Figure(; fontsize=FONTSIZE)
    ax = Axis(fig[1, 1]; xlabel=L"$\sqrt{N}$", yscale=log10, ylabel=L"$$relative error")
    ylims!(ax, 1e-16, 1)

    sas = src.vicsek_2d(1 / 3)
    k = 5.0
    Np = 750
    f_diam = 1 / 100

    colors, i = Makie.to_colormap(:tab10), 0
    for xs in xs_list
        for x in xs
            i += 1
            y = SVector(x, -0.1)
            fct = z -> green_kernel(k, z, y)

            result, precision = reference_h(fct; sas=sas, nb_pts_cbt=15, f_diam=f_diam / 2)
            @show precision

            nb_pts, h, val = sequence_p_version(fct; sas=sas, nb_pts_max=Np)

            scatterlines!(ax, sqrt.(nb_pts), relative_error.(val, result); color=colors[i])
        end

        if SAVEFIG
            save("2d-vicsek-singular-pv-$i.pdf", fig)
        end
    end

    return fig
end

# ╔═╡ 946a27ef-d769-4251-b032-85bba79f6a16
plot_vicsek_pv_sing([[-2.0], [-1.4, -1.2], [-1.0, -0.67]])

# ╔═╡ fa0b805f-6b2e-4828-bded-6c150e95b5d0
function plot_vicsek_hv_sing(xs_list::Vector{Vector{Float64}})
    fig = Figure(; fontsize=FONTSIZE)
    ax = Axis(
        fig[1, 1]; xscale=log10, xlabel=L"$h$", yscale=log10, ylabel=L"$$relative error"
    )
    ylims!(ax, 1e-16, 1)

    sas = src.vicsek_2d(1 / 3)
    k = 5.0
    f_diam = 1 / 100

    colors, i = Makie.to_colormap(:tab10), 0
    for xs in xs_list
        for x in xs
            i += 1
            y = SVector(x, -0.1)
            fct = z -> green_kernel(k, z, y)

            result, precision = reference_h(fct; sas=sas, nb_pts_cbt=15, f_diam=f_diam / 2)
            @show precision

            for deg in (3, 5)
                nb_pts, h, val = sequence_h_version(
                    fct;
                    sas=sas,
                    cbt=src.compute_cubature(sas, POINTTYPE, deg; maxiter=MAXITER),
                    f_diam,
                )
                scatterlines!(
                    ax, h, relative_error.(val, result); color=colors[i], linestyle=:dash
                )
            end
        end

        if SAVEFIG
            save("2d-vicsek-singular-hv-$i.pdf", fig)
        end
    end

    return fig
end

# ╔═╡ 4e1fac1c-74f2-48d1-8548-4de3451f150a
plot_vicsek_hv_sing([[-2.0], [-1.4, -1.2], [-1.0, -0.67]])

# ╔═╡ 84b0789e-c142-45d0-8b26-d806a602ebbc
function plot_vicsek_pv()
    fig = Figure(; fontsize=FONTSIZE)
    ax = Axis(fig[1, 1]; xlabel=L"$\sqrt{N}$", yscale=log10, ylabel=L"$$relative error")
    ylims!(ax, 1e-16, 1)

    k = 5.0
    Np = 750
    f_diam = 1 / 100

    for (sas, c) in zip(
        (src.vicsek_2d(1 / 3), src.vicsek_2d(1 / 3, 0.4), src.vicsek_2d(1 / 3, π / 4)),
        Makie.to_colormap(:tab10),
    )
        y = SVector(-2.0, -0.1)
        fct = z -> green_kernel(k, z, y)

        result, precision = reference_h(fct; sas=sas, nb_pts_cbt=15, f_diam=f_diam / 2)
        @show precision

        nb_pts, h, val = sequence_p_version(fct; sas=sas, nb_pts_max=Np)

        scatterlines!(ax, sqrt.(nb_pts), relative_error.(val, result); color=c)
    end

    if SAVEFIG
        save("2d-vicsek-pv.pdf", fig)
    end

    return fig
end

# ╔═╡ 656b1413-ba8f-4913-a397-f4271e5fdcf6
plot_vicsek_pv()

# ╔═╡ 1b418087-c6b6-458c-a4e0-02a7abe23e8b
function plot_vicsek_hv()
    fig = Figure(; fontsize=FONTSIZE)
    ax = Axis(
        fig[1, 1]; xscale=log10, xlabel=L"$h$", yscale=log10, ylabel=L"$$relative error"
    )
    ylims!(ax, 1e-16, 1)

    sas = src.vicsek_2d(1 / 3)
    k = 5.0
    f_diam = 1 / 100

    for (sas, c) in zip(
        (src.vicsek_2d(1 / 3), src.vicsek_2d(1 / 3, 0.4), src.vicsek_2d(1 / 3, π / 4)),
        Makie.to_colormap(:tab10),
    )
        y = SVector(-2.0, -0.1)
        fct = z -> green_kernel(k, z, y)

        result, precision = reference_h(fct; sas=sas, nb_pts_cbt=15, f_diam=f_diam / 2)
        @show precision

        for deg in (3, 5)
            nb_pts, h, val = sequence_h_version(
                fct;
                sas=sas,
                cbt=src.compute_cubature(sas, POINTTYPE, deg; maxiter=MAXITER),
                f_diam,
            )
            scatterlines!(ax, h, relative_error.(val, result); color=c, linestyle=:dash)
        end
    end

    if SAVEFIG
        save("2d-vicsek-hv.pdf", fig)
    end

    return fig
end

# ╔═╡ 0f6d9df4-da16-43c0-876f-a1646393ce59
plot_vicsek_hv()

# ╔═╡ Cell order:
# ╠═a6000660-39a4-11f0-139e-59028a2cdaa8
# ╠═190c44f9-f13c-4d6d-be84-b250e0f3ff61
# ╠═26308e9f-52b5-465b-b121-96ee8c50aab3
# ╠═0a1f6f87-f80c-4aaf-9162-3495b00152c1
# ╠═946a27ef-d769-4251-b032-85bba79f6a16
# ╠═fe3d2924-9fd2-4fd9-b0e3-6522e1d95edf
# ╠═4e1fac1c-74f2-48d1-8548-4de3451f150a
# ╠═fa0b805f-6b2e-4828-bded-6c150e95b5d0
# ╠═656b1413-ba8f-4913-a397-f4271e5fdcf6
# ╠═84b0789e-c142-45d0-8b26-d806a602ebbc
# ╠═0f6d9df4-da16-43c0-876f-a1646393ce59
# ╠═1b418087-c6b6-458c-a4e0-02a7abe23e8b
# ╠═0a613021-6729-484d-b919-7690a45c6493
# ╠═033139d9-6c56-4e10-924f-306c0b2fa627
# ╠═2c55c152-b2e3-432a-a781-b6045f12e5e2
# ╠═1641cc06-e7f2-4d39-89fe-0446c6edbfe9
# ╠═358da4c0-022b-4232-ad23-9e1e34729463
# ╠═21180f74-48bf-4a85-a818-3127799681d4
