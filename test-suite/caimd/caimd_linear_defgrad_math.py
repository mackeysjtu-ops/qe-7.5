#!/usr/bin/env python3
"""Deterministic Phase 1 tests for CAIMD linear_defgrad mathematics.

These tests intentionally do not run scientific AIMD.  They mirror the formulas
implemented in PW/src/caimd_module.f90 and also check that the production source
contains the required F_load, A_event, and incremental current-cell operations.
"""

from __future__ import annotations

import argparse
import math
from pathlib import Path


TOL = 1.0e-14
DET_FLOOR = 1.0e-18


def eye() -> list[list[float]]:
    return [[1.0 if i == j else 0.0 for j in range(3)] for i in range(3)]


def zeros() -> list[list[float]]:
    return [[0.0 for _ in range(3)] for _ in range(3)]


def matmul(a: list[list[float]], b: list[list[float]]) -> list[list[float]]:
    return [[sum(a[i][k] * b[k][j] for k in range(3)) for j in range(3)] for i in range(3)]


def transpose(a: list[list[float]]) -> list[list[float]]:
    return [[a[j][i] for j in range(3)] for i in range(3)]


def det3(a: list[list[float]]) -> float:
    return (
        a[0][0] * (a[1][1] * a[2][2] - a[1][2] * a[2][1])
        - a[0][1] * (a[1][0] * a[2][2] - a[1][2] * a[2][0])
        + a[0][2] * (a[1][0] * a[2][1] - a[1][1] * a[2][0])
    )


def inv3(a: list[list[float]]) -> list[list[float]]:
    det = det3(a)
    if not math.isfinite(det) or abs(det) <= DET_FLOOR:
        raise ValueError("singular matrix")
    return [
        [
            (a[1][1] * a[2][2] - a[1][2] * a[2][1]) / det,
            -(a[0][1] * a[2][2] - a[0][2] * a[2][1]) / det,
            (a[0][1] * a[1][2] - a[0][2] * a[1][1]) / det,
        ],
        [
            -(a[1][0] * a[2][2] - a[1][2] * a[2][0]) / det,
            (a[0][0] * a[2][2] - a[0][2] * a[2][0]) / det,
            -(a[0][0] * a[1][2] - a[0][2] * a[1][0]) / det,
        ],
        [
            (a[1][0] * a[2][1] - a[1][1] * a[2][0]) / det,
            -(a[0][0] * a[2][1] - a[0][1] * a[2][0]) / det,
            (a[0][0] * a[1][1] - a[0][1] * a[1][0]) / det,
        ],
    ]


def f_load(event_index: int, stride: int, d: list[list[float]]) -> list[list[float]]:
    f = eye()
    factor = float(event_index * stride)
    for i in range(3):
        for j in range(3):
            f[i][j] += factor * d[i][j]
    return f


def green_lagrange(f: list[list[float]]) -> list[list[float]]:
    c = matmul(transpose(f), f)
    e = [[0.5 * c[i][j] for j in range(3)] for i in range(3)]
    for i in range(3):
        e[i][i] -= 0.5
    return e


def a_event(event_index: int, stride: int, d: list[list[float]]) -> list[list[float]]:
    before = f_load(event_index - 1, stride, d)
    after = f_load(event_index, stride, d)
    if det3(before) <= DET_FLOOR:
        raise ValueError("F_load before event has non-positive determinant")
    if det3(after) <= DET_FLOOR:
        raise ValueError("F_load after event has non-positive determinant")
    a = matmul(after, inv3(before))
    if det3(a) <= DET_FLOOR:
        raise ValueError("A_event has non-positive determinant")
    return a


def event_steps(stride: int, count: int) -> list[int]:
    return [1 + i * stride for i in range(count)]


def validate_controls(mode: str = "legacy_strain", step_strain: list[float] | None = None,
                      step_defgrad: list[list[float]] | None = None, stride: int = 1) -> None:
    step_strain = [0.0] * 6 if step_strain is None else step_strain
    step_defgrad = zeros() if step_defgrad is None else step_defgrad
    if stride <= 0:
        raise ValueError("caimd_stride must be positive")
    if mode not in {"legacy_strain", "linear_defgrad"}:
        raise ValueError("unknown caimd_loading_mode")
    if any(not math.isfinite(x) for x in step_strain):
        raise ValueError("non-finite caimd_step_strain")
    if any(not math.isfinite(x) for row in step_defgrad for x in row):
        raise ValueError("non-finite caimd_step_defgrad")
    if mode == "legacy_strain" and any(abs(x) > 0.0 for row in step_defgrad for x in row):
        raise ValueError("caimd_step_defgrad requires linear_defgrad")
    if mode == "linear_defgrad" and any(abs(x) > 0.0 for x in step_strain):
        raise ValueError("caimd_step_strain requires legacy_strain")


def assert_close(name: str, value: float, expected: float, tol: float = TOL) -> None:
    if abs(value - expected) > tol:
        raise AssertionError(f"{name}: got {value:.17g}, expected {expected:.17g}")


def assert_matrix_close(name: str, a: list[list[float]], b: list[list[float]], tol: float = TOL) -> None:
    for i in range(3):
        for j in range(3):
            assert_close(f"{name}[{i+1},{j+1}]", a[i][j], b[i][j], tol)


def test_source_formula_presence(root: Path) -> None:
    src = (root / "PW/src/caimd_module.f90").read_text(errors="ignore")
    required = [
        "caimd_build_f_load",
        "REAL( event_index, DP ) * REAL( caimd_stride_, DP )",
        "caimd_last_a_event_ = MATMUL( caimd_last_f_after_, f_before_inv )",
        "h_new = MATMUL( caimd_last_a_event_, h_pr )",
        "caimd_green_lagrange",
        "caimd_transform_cell_history",
        "TRIM( caimd_loading_mode_ ) == 'linear_defgrad' )",
    ]
    missing = [item for item in required if item not in src]
    if missing:
        raise AssertionError(f"production formula source checks missing: {missing}")


def run_tests(root: Path) -> None:
    validate_controls()
    test_source_formula_presence(root)

    d = zeros()
    assert_matrix_close("identity before first event", f_load(0, 10, d), eye())

    d = zeros()
    d[0][2] = 2.0e-6
    f = f_load(10000, 10, d)
    assert_close("simple shear F13", f[0][2], 0.20)
    assert_close("simple shear det", det3(f), 1.0)
    e = green_lagrange(f)
    assert_close("simple shear E13", e[0][2], 0.10)
    assert_close("simple shear E31", e[2][0], 0.10)
    assert_close("simple shear E33", e[2][2], 0.02, tol=1.0e-13)

    d = zeros()
    d[2][2] = -2.0e-6
    f = f_load(7500, 10, d)
    assert_close("compression F33", f[2][2], 0.85)
    assert_close("compression det", det3(f), 0.85)
    e = green_lagrange(f)
    assert_close("compression E33", e[2][2], -0.13875)
    a = a_event(7500, 10, d)
    expected = (1.0 - 7500 * 10 * 2.0e-6) / (1.0 - 7499 * 10 * 2.0e-6)
    assert_close("incremental compression A33", a[2][2], expected)

    for sign in (1.0, -1.0):
        d = zeros()
        d[0][2] = sign * 2.0e-6
        current = eye()
        for m in range(1, 10001):
            current = matmul(a_event(m, 10, d), current)
        assert_matrix_close(f"repeated shear sign {sign}", current, f_load(10000, 10, d), tol=5.0e-13)

    for value in (2.0e-6, -2.0e-6):
        d = zeros()
        d[2][2] = value
        assert det3(f_load(100, 10, d)) > 0.0

    d = zeros()
    d[2][2] = -0.2
    try:
        a_event(1, 10, d)
    except ValueError:
        pass
    else:
        raise AssertionError("singular or inverted compression was not rejected")

    d = zeros()
    d[0][2] = 2.0e-6
    try:
        validate_controls("legacy_strain", step_defgrad=d)
    except ValueError:
        pass
    else:
        raise AssertionError("legacy mode accepted caimd_step_defgrad")
    try:
        validate_controls("linear_defgrad", step_strain=[1.0e-6, 0.0, 0.0, 0.0, 0.0, 0.0])
    except ValueError:
        pass
    else:
        raise AssertionError("linear_defgrad accepted caimd_step_strain")
    try:
        validate_controls("not_a_mode")
    except ValueError:
        pass
    else:
        raise AssertionError("invalid loading mode was not rejected")

    h_current = [[1.01, 0.02, 0.0], [0.00, 0.98, 0.03], [0.01, 0.00, 1.04]]
    d = zeros()
    d[0][2] = 2.0e-6
    incremental = matmul(a_event(100, 10, d), h_current)
    reset_like = f_load(100, 10, d)
    if max(abs(incremental[i][j] - reset_like[i][j]) for i in range(3) for j in range(3)) < 1.0e-6:
        raise AssertionError("incremental current-cell update collapsed to reference reset")

    if event_steps(10, 4) != [1, 11, 21, 31]:
        raise AssertionError("stride-10 event convention changed")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--source-root", type=Path, default=Path.cwd())
    args = parser.parse_args()
    run_tests(args.source_root.resolve())
    print("CAIMD linear_defgrad deterministic math tests passed")


if __name__ == "__main__":
    main()
