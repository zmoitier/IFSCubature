using Printf, LinearAlgebra, StaticArrays, DataStructures, TOML, CairoMakie

import IFSCubature as src

const POINTTYPE = "Chebyshev-1"
const MAXITER = 2048

function laplace_kernel(x, y)
    return 1 / norm(x - y)
end

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

function relative_error(result::Number, reference::Number)
    if isapprox(reference, 0)
        return abs(result)
    end
    return abs(result / reference - 1)
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

k = 0.0
x0 = SVector{2,Float64}([0.1, -2.0])
if k ≈ 0
    fct = x -> laplace_kernel(x, x0)
else
    fct = x -> green_kernel(k, x, x0)
end

fct = x -> x[1] .^ 30 * x[2] .^ 30

# sas = src.cantor_dust(1 / 3, [-1.0, 1.0], 2)
# sas = cantor_dust_test()
sas = src.cantor_dust_non_sym()

fig = Figure()
ax = Axis(fig[1, 1]; yscale=log10)

v_ref, prec = reference_h(fct; sas=sas, nb_pts_cbt=15, f_diam=1 / 200)
println((v_ref, prec))
n, h, v = sequence_p_version(fct; sas=sas, nb_pts_max=1000)

lines!(n, 0.1 .* exp.((-0.25) .* n .^ 0.5); color=:black)

e = abs.(v ./ v_ref .- 1) .+ 1e-17
scatterlines!(ax, n, e; linestyle=:dash)

display(fig)
