"""Affine maps."""

from dataclasses import dataclass
from functools import reduce

from IPython.display import display
from sympy import Eq, Expr, Lt, Matrix, Ne, cos, diag, eye, sin, sqrt, zeros


@dataclass(frozen=True, slots=True)
class AffineMap:
    """Affine maps defined by x -> Ax+b."""

    A: Matrix  # (n, n) square matrix
    b: Matrix  # (n, 1) matrix aka a vector

    def __post_init__(self) -> None:
        """Check the size."""

        A_shape = self.A.shape
        if A_shape[0] != A_shape[1]:
            raise ValueError(f"A {A_shape} must be a square matrix.")

        if self.b.shape != (A_shape[1], 1):
            raise ValueError(f"b {self.b.shape} must be a ({A_shape[1]},1) matrix.")

    def __call__(self, x: Matrix) -> Matrix:
        return self.A @ x + self.b

    @property
    def space_dim(self) -> int:
        """Return the space dimension."""
        return self.A.shape[0]


def _compose(f: AffineMap, g: AffineMap) -> AffineMap:
    """Return the affine transformation of g \\circ f."""

    return AffineMap(g.A * f.A, g.A * f.b + g.b)


def compose(*affine_maps: AffineMap) -> AffineMap:
    """For a list of affine maps [f1, f2, ..., fn]
    return fn \\circ ... \\circ f2 \\circ f1."""

    return reduce(_compose, affine_maps)


class AffineContraction(AffineMap):
    """Affine contraction."""

    def check_contraction(self) -> None:
        """Check if self is a contraction."""
        for s in self.A.singular_values():
            display(Lt(s, 1))

    def fix_point(self) -> Matrix:
        """Return the fix point of the affine contraction."""
        return (eye(self.space_dim) - self.A).solve(self.b)


def contractive_similarity(
    contraction_factor: Expr,
    orthogonal_matrix: Matrix,
    fix_point: Matrix,
) -> AffineContraction:
    """Return the contractive similarity x -> rT(x-c)+c from the contraction factor `r`,
    the orthogonal matrix `T` and the fix point `c`."""

    A = contraction_factor * orthogonal_matrix
    b = fix_point - A * fix_point

    return AffineContraction(A, b)


class AffineTransformation(AffineMap):
    """Affine transformation."""

    def check_invertible(self) -> None:
        """Check if self is invertible."""
        display(Ne(self.A.det(), 0))


def translation(vector: Matrix) -> AffineTransformation:
    """Return the translation."""
    dim, m = vector.shape

    if m != 1:
        raise ValueError("vector must be an (n, 1) matrix.")

    return AffineTransformation(eye(dim), vector)


def homothety(ratio: Expr, dim: int) -> AffineTransformation:
    """Homothety of ratio = `ratio` in dimension `dim`."""

    return AffineTransformation(diag(*[ratio for _ in range(dim)]), zeros(dim, 1))


def reflection(normal: Matrix) -> AffineTransformation:
    """Reflection with respect to the hyperplane define by its normal."""

    dim, m = normal.shape

    if m != 1:
        raise ValueError("vector must be an (n, 1) matrix.")

    trans = normal.T
    norm2 = (trans @ normal)[0, 0]
    if Eq(norm2, 0):
        raise ValueError(f"{normal} must be a non-zero vector.")

    return AffineTransformation(eye(dim) - 2 * normal * trans / norm2, zeros(dim, 1))


def rotation_2d(angle: Expr) -> AffineTransformation:
    """2d Rotation of angle = `angle`."""

    return AffineTransformation(rotation_matrix_2d(angle), zeros(2, 1))


def rotation_matrix_2d(angle: Expr) -> Matrix:
    """Return the 2d rotation matrix"""

    c, s = cos(angle).trigsimp(), sin(angle).trigsimp()
    return Matrix([[c, -s], [s, c]])


def rotation_3d(angle: Expr, axis: tuple[Expr, Expr, Expr]) -> AffineTransformation:
    """3d Rotation of angle = `angle` and axis = `axis`."""

    return AffineTransformation(rotation_matrix_3d(angle, axis), zeros(3, 1))


def rotation_matrix_3d(angle: Expr, axis: tuple[Expr, Expr, Expr]) -> Matrix:
    """Return the 3d rotation matrix"""

    axis_norm = sqrt(sum(v**2 for v in axis))

    if Eq(axis_norm, 0):
        raise ValueError(f"{axis} must be a non-zero vector.")

    c, s = cos(angle).trigsimp(), sin(angle).trigsimp()

    u = Matrix(axis) / axis_norm

    cross_product_matrix = Matrix(
        [
            [0, -u[2], u[1]],
            [u[2], 0, -u[0]],
            [-u[1], u[0], 0],
        ]
    )

    outer_product = u @ u.T

    return c * eye(3) + s * cross_product_matrix + (1 - c) * outer_product
