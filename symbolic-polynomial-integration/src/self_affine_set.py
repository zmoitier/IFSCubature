"""Self-Affine set from IFS."""

from dataclasses import dataclass

from sympy import Expr, Rational

from .affine_map import AffineContraction


@dataclass(frozen=True, slots=True)
class SelfAffineSet:
    """Attractor from an IFS."""

    ifs: list[AffineContraction]
    measure: list[Expr]
    space_dim: int

    def __post_init__(self) -> None:
        if len(self.ifs) != len(self.measure):
            raise ValueError("ifs and measure must have the same length.")

        for S in self.ifs:
            if S.space_dim != self.space_dim:
                raise ValueError(f"{S} must be an {self.space_dim}d affine map.")

    @property
    def nb_maps(self) -> int:
        """Return number of contractive maps."""
        return len(self.ifs)

    @property
    def space_dim(self) -> int:
        """Return the space dimension."""
        return self.ifs[0].space_dim


def uniform_measure(n: int) -> list[Expr]:
    """Return the n-tuple (1/n, ..., 1/n)."""

    return [Rational(1, n) for _ in range(n)]
