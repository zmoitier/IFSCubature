"""Export values to file."""

from pathlib import Path
from typing import Any

from sympy import Expr

from .self_affine_set import SelfAffineSet


def _dumps_value(value: Any) -> str:
    """Recursively dump value."""

    match value:
        case bool():
            return "true" if value else "false"

        case float() | int():
            return str(value)

        case str():
            return f'"{value}"'

        case Expr():
            return str(value.evalf(17))

        case list() | tuple():
            return f"[{', '.join(_dumps_value(v) for v in value)}]"

        case dict():
            key_val = ", ".join(
                f'"{key}" = {_dumps_value(val)}' for key, val in value.items()
            )
            return f"{{ {key_val} }}"

        case _:
            raise TypeError(f"{type(value).__name__} is not supported")


def _dumps(toml_dict: dict[str, Any], table: str) -> str:
    """Dump in str."""

    document = []
    for key, value in toml_dict.items():
        match value:
            case dict():
                table_key = f"{table}.{key}" if table else key
                document.append(f"\n[{table_key}]")
                document.append(_dumps(value, table_key))

            case _:
                document.append(f"{key} = {_dumps_value(value)}")

    return "\n".join(document)


def export_to_file(
    values: dict[tuple[int, ...], Expr], attractor: SelfAffineSet, filename: Path | str
) -> None:
    """Export polynomials integral to file."""

    toml_dict: dict[str, Any] = {}

    toml_dict["space-dim"] = attractor.space_dim
    toml_dict["ifs"] = [{"A": S.A.tolist(), "b": S.b.tolist()} for S in attractor.ifs]
    toml_dict["measure"] = attractor.measure

    toml_dict["exponent-to-value"] = {
        "-".join(str(e) for e in key): val for key, val in values.items()
    }

    toml = _dumps(toml_dict, "")

    with open(f"{filename}.toml", "w") as file:
        file.write(toml)

    return None
