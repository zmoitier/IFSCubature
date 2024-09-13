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
function plot_eigenvalues(sas::src.SelfAffineSet{D,T,N},nb_pts_cbt::Int) where {D,T,N}
S, _ = src.matrix_from_lagrange_polynomials(sas, POINTTYPE, nb_pts_cbt)
	eigs = eigvals(S)

	fig = Figure(; fontsize=FONTSIZE)
	ax = PolarAxis(fig[1,1])

    scatter!(ax, abs.(eigs), eigs)

    return fig

end

# ╔═╡ 51bf816b-d891-4442-bfe6-ff095baa950b
plot_eigenvalues(src.cantor_set(1/3,[0.0,1.0]), 32)

# ╔═╡ Cell order:
# ╠═81f51946-7054-11ef-3aa5-29727a98a16c
# ╠═458c751e-54d3-4d28-a713-d779f7d167e4
# ╠═51bf816b-d891-4442-bfe6-ff095baa950b
# ╠═817364c0-24fa-4827-b2ab-ba301db71f67
