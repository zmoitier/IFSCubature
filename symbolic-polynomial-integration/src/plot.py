"""Plot utilities."""

from matplotlib.pyplot import subplots
from numpy import arange, array, float64, zeros, zeros_like
from numpy.random import choice
from numpy.typing import NDArray
from sympy import matrix2numpy

from .self_affine_set import SelfAffineSet


def render_chaos_game(attractor: SelfAffineSet, nb_points: int) -> None:
    """Render the attractor using chaos game method."""

    fix_points = [matrix2numpy(S.fix_point(), dtype=float64) for S in attractor.ifs]

    ifs = [
        (matrix2numpy(S.A, dtype=float64), matrix2numpy(S.b, dtype=float64))
        for S in attractor.ifs
    ]

    probability = array([p for p in attractor.measure], dtype=float64)

    n = attractor.space_dim

    range_L = arange(attractor.nb_maps)

    M = nb_points // attractor.nb_maps

    if n == 3:
        _, ax = subplots(subplot_kw={"projection": "3d"}, layout="constrained")
    else:
        _, ax = subplots(layout="constrained")

    for c in fix_points:
        coords = zeros((attractor.space_dim, M), dtype=float64)
        v = c

        coords[:, 0] = v[:, 0]
        for m in range(1, M):
            i = choice(range_L, p=probability)
            A, b = ifs[i]
            v = A @ v + b
            coords[:, m] = v[:, 0]

        match n:
            case 1:
                ax.plot(coords[0, :], zeros_like(coords[0, :]), "C0,")
            case 2:
                ax.plot(coords[0, :], coords[1, :], "C0,")
            case 3:
                ax.plot(coords[0, :], coords[1, :], coords[2, :], "C0,")
            case _:
                raise ValueError(f"no visual for spacial dimension {n}.")

    ax.set_aspect("equal", "box")
    ax.grid(True)


def render_pre_attractor(attractor: SelfAffineSet, start: NDArray, level: int) -> None:
    """Render pre-attractor up to level `level` from `start`."""

    ifs = [
        (matrix2numpy(S.A, dtype=float64), matrix2numpy(S.b, dtype=float64))
        for S in attractor.ifs
    ]

    match n := attractor.space_dim:
        case 1:
            _render_pre_1d(ifs, start, level)
        case 2:
            _render_pre_2d(ifs, start, level)
        case 3:
            _render_pre_3d(ifs, start, level)
        case _:
            raise ValueError(f"no visual for spacial dimension {n}.")

    return None


def _render_pre_1d(
    ifs: list[tuple[NDArray, NDArray]], start: NDArray, level: int
) -> None:
    """Render 1d pre-attractor up to level `level` from `start`."""

    _, axs = subplots(nrows=level + 1, figsize=[6.4, 4.8], layout="constrained")

    shapes = [start]

    for segment in shapes:
        axs[0].plot(segment[0, :], zeros(2), "C0-", lw=2)

    axs[0].set_aspect("equal", "box")
    axs[0].grid(True)

    for p in range(1, level + 1):
        shapes = _apply_ifs(ifs, shapes)

        for segment in shapes:
            axs[p].plot(segment[0, :], zeros(2), "C0-", lw=2)

        axs[p].set_aspect("equal", "box")
        axs[p].grid(True)

    return None


def _render_pre_2d(
    ifs: list[tuple[NDArray, NDArray]], start: NDArray, level: int
) -> None:
    """Render 2d pre-attractor up to level `level` from `start`."""

    _, axs = subplots(ncols=level + 1, figsize=[level * 6.4, 4.8], layout="constrained")

    shapes = [start]
    kwargs = {
        "facecolor": (0.5, 0.5, 0.5, 0.5),
        "edgecolor": "k",
        "linewidth": 2,
        "zorder": 2,
    }

    for shape in shapes:
        axs[0].fill(shape[0, :], shape[1, :], **kwargs)

    axs[0].set_aspect("equal", "box")
    axs[0].grid(True, zorder=1)

    for p in range(1, level + 1):
        shapes = _apply_ifs(ifs, shapes)

        for shape in shapes:
            axs[p].fill(shape[0, :], shape[1, :], **kwargs)

        axs[p].set_aspect("equal", "box")
        axs[p].grid(True, zorder=1)

    return None


def _render_pre_3d(
    ifs: list[tuple[NDArray, NDArray]], start: NDArray, level: int
) -> None:
    """Render 3d pre-attractor up to level `level` from `start`."""

    _, axs = subplots(
        ncols=level + 1, subplot_kw={"projection": "3d"}, layout="constrained"
    )

    shapes = [start]

    for segment in shapes:
        axs[0].plot(segment[0, :], zeros(2), "C0-", lw=2)

    axs[0].set_aspect("equal", "box")
    axs[0].grid(True)

    for p in range(1, level + 1):
        shapes = _apply_ifs(ifs, shapes)

        for segment in shapes:
            axs[p].plot(segment[0, :], zeros(2), "C0-", lw=2)

        axs[p].set_aspect("equal", "box")
        axs[p].grid(True)

    return None


def _apply_ifs(
    ifs: list[tuple[NDArray, NDArray]], shapes: list[NDArray]
) -> list[NDArray]:
    """Apply IFS to pre-attractor."""

    tmp = []
    for shape in shapes:
        for fct in ifs:
            A, b = fct
            tmp.append(A @ shape + b)

    return tmp
