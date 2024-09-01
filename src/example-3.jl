"""Return the cube [a, b]x[a, b]x[a, b]."""
function cube(a::T, b::T) where {T}
    return cantor_dust(T(1//2), [a, b], 3, "3d-cube")
end

# """Return the tetrahedron."""
# function tetrahedron()
# end

"""Return the Sierpinski tetrahedron."""
function sierpinski_tetrahedron(contraction_factor::Union{Float64,Nothing}=nothing)
    if isnothing(contraction_factor)
        ρ = 0.5
    else
        ρ = contraction_factor
    end

    @assert 0 ≤ ρ ≤ 1 "the contraction factor must be between [0, 1)."

    ifs = [
        contractive_similarity(ρ, c) for
        c in [[1, 0, -√2 / 2], [-1, 0, -√2 / 2], [0, 1, √2 / 2], [0, -1, √2 / 2]]
    ]

    b_ball = hyper_ball(fill(0.0, 3), √(3 / 2))
    b_box = hyper_box(fill(0.0, 3), Matrix(Diagonal([1.0, 1.0, √2 / 2])))

    return SelfAffineSet(ifs, fill(1 / 4, 4), b_ball, b_box, "3d-sierpinski-tetrahedron")
end

"""Return Menger sponge."""
function menger_sponge()
    ρ = 1 / 3
    fix_points = [
        [1, 0, 1],
        [1, 1, 1],
        [0, 1, 1],
        [-1, 1, 1],
        [-1, 0, 1],
        [-1, -1, 1],
        [0, -1, 1],
        [1, -1, 1],
        #
        [1, 1, 0],
        [-1, 1, 0],
        [-1, -1, 0],
        [1, -1, 0],
        #
        [1, 0, -1],
        [1, 1, -1],
        [0, 1, -1],
        [-1, 1, -1],
        [-1, 0, -1],
        [-1, -1, -1],
        [0, -1, -1],
        [1, -1, -1],
    ]

    ifs = [contractive_similarity(ρ, c) for c in fix_points]

    b_ball = hyper_ball(fill(0.0, 3), √3)
    b_box = hyper_box(fill(0.0, 3), Matrix(Diagonal(fill(1.0, 3))))

    return SelfAffineSet(ifs, fill(1 / 20, 20), b_ball, b_box, "3d-menger-sponge")
end

"""Return the 3d Vicsek fractal."""
function vicsek_3d(contraction_factor::Float64, central_rotation::Bool=false)
    @assert 0 ≤ contraction_factor ≤ 1 "the contraction factor must be between [0, 1)."

    name = "3d-vicsek"
    if central_rotation
        Ry = matrix_rotation_3d(atan(√2), [0.0, 1.0, 0.0])
        Rz = matrix_rotation_3d(-0.25, [0.0, 0.0, 1.0]; implicit_pi=true)
        R = Rz * Ry
        name *= "-rot"
    else
        R = Matrix(I, 3, 3)
    end
    S0 = contractive_similarity(contraction_factor, R, fill(0.0, 3))

    ifs::Vector{AffineMap{3,Float64,9}} = [S0]
    for c in product([-1, 1], [-1, 1], [-1, 1])
        push!(ifs, contractive_similarity(contraction_factor, [v for v in c]))
    end

    b_ball = hyper_ball(fill(0.0, 3), √3)
    b_box = hyper_box(fill(0.0, 3), Matrix(Diagonal(fill(1.0, 3))))

    return SelfAffineSet(ifs, fill(1 / 9, 9), b_ball, b_box, name)
end
