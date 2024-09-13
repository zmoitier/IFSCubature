### A Pluto.jl notebook ###
# v0.19.46

using Markdown
using InteractiveUtils

# ╔═╡ 81f51946-7054-11ef-3aa5-29727a98a16c
begin
    using Pkg: Pkg
    Pkg.activate(Base.current_project())
    Pkg.instantiate()

    using LinearAlgebra, CairoMakie

    import IFSCubature as src
end

# ╔═╡ 458c751e-54d3-4d28-a713-d779f7d167e4
begin
    #! Type of points
    # const POINTTYPE = "equispaced"
    const POINTTYPE = "Chebyshev-1"
    # const POINTTYPE = "Chebyshev-2"
    # const POINTTYPE = "Gauss-Legendre"
    # const POINTTYPE = "Gauss-Lobatto"

    const MAXITER = 2048

    #! Ploting constants
    const FONTSIZE = 20

    "Global parameters"
end

# ╔═╡ 817364c0-24fa-4827-b2ab-ba301db71f67
function plot_eigenvalues(sas::src.SelfAffineSet{D,T,N}, nb_pts_cbt::Int) where {D,T,N}
    S, _ = src.matrix_from_lagrange_polynomials(sas, POINTTYPE, nb_pts_cbt)
    eigs = eigvals(S)

    fig = Figure(; size=(800, 450), fontsize=FONTSIZE)
    axs = [
        PolarAxis(fig[1, 1]; rlimits=(0, 1.05), rticks=[0, 0.5, 1], title=sas.name),
        PolarAxis(fig[1, 2]; title=L"(\log(r),\ \theta)"),
    ]

    n = sum(isapprox.(eigs, 1))
    println("number of eigenvalues ≈ 1 = $n")

    r = abs.(eigs)
    t = angle.(eigs)

    scatter!(axs[1], t, r; marker=:cross, markersize=12)
    scatter!(axs[2], t, log10.(r); marker=:cross, markersize=12)

    return fig
end

# ╔═╡ 51bf816b-d891-4442-bfe6-ff095baa950b
plot_eigenvalues(src.cantor_set(0.5, [0.0, 1.0]), 32)

# ╔═╡ f2ba5941-ea44-4925-bbf7-50c22ff2d0ad
begin
    ρ = [0.1, 0.2, 0.3]
    d = src.similarity_dimension(ρ)
    sas = src.cantor_set(ρ, [-1.0, 0.0, 1.0], ρ .^ d)

    plot_eigenvalues(sas, 32)
end

# ╔═╡ d18583bb-5772-449b-ae16-20d2816598bb
plot_eigenvalues(src.cantor_dust(1 / 3, [0.0, 1.0], 2), 16)

# ╔═╡ 578405ae-82d0-46f1-bc20-6ddb3a2da1f2
plot_eigenvalues(src.sierpinski_triangle(), 16)

# ╔═╡ f5845ab8-99be-4396-af2d-b151a9fe5c24
plot_eigenvalues(src.sierpinski_triangle_fat(2), 16)

# ╔═╡ bbef60c3-584b-4c6f-be12-1bd0625f7457
plot_eigenvalues(src.vicsek_2d(1 / 3), 16)

# ╔═╡ 1cf6eb95-aa34-4433-b771-89b1f27bd850
plot_eigenvalues(src.vicsek_2d(1 / 3, 0.4), 16)

# ╔═╡ e8f7ed8d-18e0-476b-8c20-8198d1083d75
plot_eigenvalues(src.vicsek_2d(1 / 3, π / 4), 16)

# ╔═╡ 2bf71314-dc21-4ce1-9e2b-0c5f6486867a
plot_eigenvalues(src.sierpinski_carpet(), 16)

# ╔═╡ c69afd93-7f66-4bcc-a4f8-8cc8a15c3918
plot_eigenvalues(src.koch_snowflake(), 16)

# ╔═╡ 310ae65d-1a52-49d1-b79f-bf1f1d060f5d
plot_eigenvalues(src.gosper_flowsnake(), 16)

# ╔═╡ 7a91b178-732e-464c-9148-fea04cc18f60
plot_eigenvalues(src.brick_2d(), 16)

# ╔═╡ 3196d771-4395-4650-9c4e-a96de0454c6c
plot_eigenvalues(src.durer_pentagon(), 16)

# ╔═╡ ac772316-e4b5-427d-8391-7ff6b276146f
plot_eigenvalues(src.fudgeflake(), 16)

# ╔═╡ ada05b50-f4ac-4633-9e25-7177a7a3d0e0
plot_eigenvalues(src.heighway_dragon(), 16)

# ╔═╡ 39ad96ab-1b27-471e-8117-c9826a790b46
plot_eigenvalues(src.levy_dragon(), 16)

# ╔═╡ 41b52d22-4222-4e70-b7b2-bcb0e50fbb2c
plot_eigenvalues(src.terdragon(), 16)

# ╔═╡ 0c5ce9ef-0ade-41cc-a1a6-422eaef1c479
plot_eigenvalues(src.twindragon(), 16)

# ╔═╡ dfe60e90-f05a-4386-8680-fa378a73ce8b
plot_eigenvalues(src.cantor_dust_non_sym(), 16)

# ╔═╡ 071c2d88-b677-4311-8e81-b20aaf6e619c
plot_eigenvalues(src.barnsley_fern(), 16)

# ╔═╡ 1e26797d-5dc5-4afb-8712-5b18eef91da7
plot_eigenvalues(src.sierpinski_tetrahedron(), 8)

# ╔═╡ 6d609263-89f6-4a58-be24-04b779412263
plot_eigenvalues(src.cantor_dust(1 / 3, [0.0, 1.0], 3), 8)

# ╔═╡ da7585f3-0226-40b0-82ad-972408e742c9
plot_eigenvalues(src.menger_sponge(), 8)

# ╔═╡ 67b90f32-24a3-457c-bd50-f510679839fb
plot_eigenvalues(src.vicsek_3d(1 / 3), 8)

# ╔═╡ fb2be838-b4fd-4d92-bf9e-98001b687f11
plot_eigenvalues(src.vicsek_3d(1 / 3, true), 8)

# ╔═╡ 423c6653-5185-4a96-b9e3-2ef2c0ed0235
plot_eigenvalues(src.brick_3d(), 8)

# ╔═╡ Cell order:
# ╠═81f51946-7054-11ef-3aa5-29727a98a16c
# ╠═458c751e-54d3-4d28-a713-d779f7d167e4
# ╠═51bf816b-d891-4442-bfe6-ff095baa950b
# ╠═f2ba5941-ea44-4925-bbf7-50c22ff2d0ad
# ╠═d18583bb-5772-449b-ae16-20d2816598bb
# ╠═578405ae-82d0-46f1-bc20-6ddb3a2da1f2
# ╠═f5845ab8-99be-4396-af2d-b151a9fe5c24
# ╠═bbef60c3-584b-4c6f-be12-1bd0625f7457
# ╠═1cf6eb95-aa34-4433-b771-89b1f27bd850
# ╠═e8f7ed8d-18e0-476b-8c20-8198d1083d75
# ╠═2bf71314-dc21-4ce1-9e2b-0c5f6486867a
# ╠═c69afd93-7f66-4bcc-a4f8-8cc8a15c3918
# ╠═310ae65d-1a52-49d1-b79f-bf1f1d060f5d
# ╠═7a91b178-732e-464c-9148-fea04cc18f60
# ╠═3196d771-4395-4650-9c4e-a96de0454c6c
# ╠═ac772316-e4b5-427d-8391-7ff6b276146f
# ╠═ada05b50-f4ac-4633-9e25-7177a7a3d0e0
# ╠═39ad96ab-1b27-471e-8117-c9826a790b46
# ╠═41b52d22-4222-4e70-b7b2-bcb0e50fbb2c
# ╠═0c5ce9ef-0ade-41cc-a1a6-422eaef1c479
# ╠═dfe60e90-f05a-4386-8680-fa378a73ce8b
# ╠═071c2d88-b677-4311-8e81-b20aaf6e619c
# ╠═1e26797d-5dc5-4afb-8712-5b18eef91da7
# ╠═6d609263-89f6-4a58-be24-04b779412263
# ╠═da7585f3-0226-40b0-82ad-972408e742c9
# ╠═67b90f32-24a3-457c-bd50-f510679839fb
# ╠═fb2be838-b4fd-4d92-bf9e-98001b687f11
# ╠═423c6653-5185-4a96-b9e3-2ef2c0ed0235
# ╠═817364c0-24fa-4827-b2ab-ba301db71f67
