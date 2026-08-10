#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
PYTHON_BIN="${PYTHON_BIN:-/Users/mengqicheng/.cache/codex-runtimes/codex-primary-runtime/dependencies/python/bin/python3}"
PW_BIN="${PW_BIN:-/private/tmp/qe_diamond8_af7ebfc_build/bin/pw.x}"

"${PYTHON_BIN}" -B "${ROOT}/scripts/four_dof_workflow.py" --pw "${PW_BIN}" prepare
"${PYTHON_BIN}" -B "${ROOT}/scripts/four_dof_workflow.py" --pw "${PW_BIN}" run-equil
"${PYTHON_BIN}" -B "${ROOT}/scripts/four_dof_workflow.py" --pw "${PW_BIN}" run-matrix
"${PYTHON_BIN}" -B "${ROOT}/scripts/four_dof_workflow.py" --pw "${PW_BIN}" analyze
"${PYTHON_BIN}" -B "${ROOT}/scripts/four_dof_workflow.py" --pw "${PW_BIN}" run-regression
