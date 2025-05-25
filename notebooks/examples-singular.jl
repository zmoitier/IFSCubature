### A Pluto.jl notebook ###
# v0.20.8

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
    # const POINTTYPE = "equispaced"
    const POINTTYPE = "Chebyshev-1"
    # const POINTTYPE = "Chebyshev-2"
    # const POINTTYPE = "Gauss-Legendre"
    # const POINTTYPE = "Gauss-Lobatto"

    const MAXITER = 2000

    #! Ploting constants
    const FONTSIZE = 20
    const SAVEFIG = false

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

# ╔═╡ 793d20ce-9b72-42f6-8ed1-b240fdfb918e
function plot_vicsek(x::Real, pts_chaos_nb::Int, α_attractor::Real, suffix::String="")
    fig = Figure(; size=(600, 400), fontsize=FONTSIZE)
    ax = Axis(fig[1, 1]; aspect=DataAspect(), xlabel=L"x", ylabel=L"y")

    sas = src.vicsek_2d(1 / 3)
    plot_chaos_game!(ax, sas, pts_chaos_nb, α_attractor)

    scatter!(ax, [x], [0]; color=:black, markersize=8)

    limits!(ax, -2.05, 1.05, -1.05, 1.05)

    if SAVEFIG && !isempty(name)
        save("2d-vicsek-sing-$suffix.pdf", fig)
    end

    return fig
end

# ╔═╡ 26308e9f-52b5-465b-b121-96ee8c50aab3
plot_vicsek(-2, 200_000, 0.25, "2.0")

# ╔═╡ a590c1e8-127c-4be4-859a-692ce2fe383d
plot_vicsek(-1.4, 200_000, 0.25, "1.4")

# ╔═╡ da024b91-046a-41ff-b712-47545ef4405b
plot_vicsek(-1.2, 200_000, 0.25, "1.2")

# ╔═╡ 4a1d8034-3e49-4ac1-9b18-9f8daec022f3
plot_vicsek(-1, 200_000, 0.25, "1.0")

# ╔═╡ 28a5aee7-3475-4bb8-adf6-d045e703aa14
plot_vicsek(-0.7, 200_000, 0.25, "0.7")

# ╔═╡ 033139d9-6c56-4e10-924f-306c0b2fa627
function green_kernel(k, x, y)
    nxy = norm(x - y)
    return exp(im * k * nxy) / nxy
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

# ╔═╡ 289fa304-d3de-4283-8fda-cbc9e49296a8
function plot_vicsek_pv_sing()
    fig = Figure(; fontsize=FONTSIZE)
    ax = Axis(fig[1, 1]; yscale=log10)

    sas = src.vicsek_2d(1 / 3)
    k = 5.0
    Np = 750
    f_diam = 1 / 100

    for x in [-2.0, -1.4, -1.2, -1, -0.7]
        y = SVector(x, 0.0)
        fct = z -> green_kernel(k, z, y)

        result, precision = reference_h(fct; sas=sas, nb_pts_cbt=15, f_diam=f_diam / 2)

        nb_pts, h, val = sequence_p_version(fct; sas=sas, nb_pts_max=Np)

        scatterlines!(ax, nb_pts, relative_error.(val, result))
    end

    return fig
end

# ╔═╡ 946a27ef-d769-4251-b032-85bba79f6a16
plot_vicsek_pv_sing()

# ╔═╡ Cell order:
# ╠═a6000660-39a4-11f0-139e-59028a2cdaa8
# ╠═190c44f9-f13c-4d6d-be84-b250e0f3ff61
# ╠═26308e9f-52b5-465b-b121-96ee8c50aab3
# ╠═a590c1e8-127c-4be4-859a-692ce2fe383d
# ╠═da024b91-046a-41ff-b712-47545ef4405b
# ╠═4a1d8034-3e49-4ac1-9b18-9f8daec022f3
# ╠═28a5aee7-3475-4bb8-adf6-d045e703aa14
# ╠═946a27ef-d769-4251-b032-85bba79f6a16
# ╠═793d20ce-9b72-42f6-8ed1-b240fdfb918e
# ╠═289fa304-d3de-4283-8fda-cbc9e49296a8
# ╠═0a613021-6729-484d-b919-7690a45c6493
# ╠═033139d9-6c56-4e10-924f-306c0b2fa627
# ╠═2c55c152-b2e3-432a-a781-b6045f12e5e2
# ╠═1641cc06-e7f2-4d39-89fe-0446c6edbfe9
# ╠═358da4c0-022b-4232-ad23-9e1e34729463
# ╠═21180f74-48bf-4a85-a818-3127799681d4
