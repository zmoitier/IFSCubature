"""Polynomial integration."""

from IPython.display import Latex, display
from sympy import Eq, Expr, Matrix, Poly, S, latex, linsolve, prod, symbols

from .self_affine_set import SelfAffineSet


def apply_op(p: Poly, attractor: SelfAffineSet) -> Poly:
    """Apply operator F."""

    X = Matrix(p.gens)
    exponent_to_coef = p.as_dict()

    q = S(0)
    for S_l, mu_l in zip(attractor.ifs, attractor.measure, strict=True):
        Y = S_l(X)
        for exponent, coef in exponent_to_coef.items():
            q += (
                mu_l * coef * prod(v**e for v, e in zip(Y, exponent, strict=True))
            ).expand()

    return Poly(q, p.gens)


def compute_polynomial_integral(
    attractor: SelfAffineSet, tot_deg_max: int
) -> dict[tuple[int, ...], Expr]:
    """Compute polynomial integral IFS attractor."""

    if tot_deg_max < 0:
        raise ValueError(f"total_degree_max={tot_deg_max} must be non-negative.")

    n = attractor.space_dim
    exponent_to_value: dict[tuple[int, ...], Expr] = {tuple(0 for _ in range(n)): S(1)}

    if tot_deg_max == 0:
        return exponent_to_value

    monomial_to_value: dict[Expr, Expr] = {}

    x = symbols(f"x0:{n}")
    X = Matrix(x)

    L = len(attractor.ifs)
    SX = [S.A @ X + S.b for S in attractor.ifs]

    for td in range(1, tot_deg_max + 1):
        exponents = _tot_deg(td, n)

        monomials = [
            prod(v**e for v, e in zip(x, exponent, strict=True))
            for exponent in exponents
        ]
        variables = [
            symbols(f"v_{''.join(str(e) for e in exponent)}") for exponent in exponents
        ]

        for monomial, variable in zip(monomials, variables, strict=True):
            monomial_to_value[monomial] = variable

        system = []
        for variable, exponent in zip(variables, exponents, strict=True):
            rhs = variable
            lhs = sum(
                (
                    (
                        attractor.measure[ell]
                        * prod(v**e for v, e in zip(SX[ell], exponent, strict=True))
                    ).expand()
                    for ell in range(L)
                ),
                start=S(0),
            )

            system.append(Eq(rhs, lhs.subs(reversed(list(monomial_to_value.items())))))

        (solution,) = linsolve(system, variables)

        for exponent, monomial, expr in zip(
            exponents, monomials, solution, strict=True
        ):
            exponent_to_value[exponent] = expr.factor()
            monomial_to_value[monomial] = expr

    return exponent_to_value


def _tot_deg(tot_deg: int, space_dim: int) -> list[tuple[int, ...]]:
    """Return the list of total degree tot_deg in space_dim."""

    if space_dim <= 0:
        return []

    if space_dim == 1:
        return [(tot_deg,)]

    return [
        (d, *t)
        for d in range(tot_deg, -1, -1)
        for t in _tot_deg(tot_deg - d, space_dim - 1)
    ]


def display_values(
    values: dict[tuple[int, ...], Expr],
    space_dim: int,
    *,
    variables: tuple[Expr, ...] | None = None,
    expand: bool = False,
) -> None:
    """Display values."""

    if variables is None:
        variables = tuple(symbols(f"x_1:{space_dim+1}", real=True))

    if space_dim == 1:
        measure = r"\mathrm{d}\mu(x)"
    else:
        measure = r"\mathrm{d}\mu(\boldsymbol{x})"

    for exponent, value in values.items():
        integrand = latex(prod(x**e for x, e in zip(variables, exponent, strict=True)))
        lhs = rf"\int_\Gamma {integrand}\, {measure}"

        rhs = value.expand() if expand else value

        display(Latex(rf"$\displaystyle {lhs} = {latex(rhs)} $"))
