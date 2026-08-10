#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")/.."
python3 scripts/controls_abc_workflow.py prepare
python3 scripts/controls_abc_workflow.py run-all
python3 scripts/controls_abc_workflow.py analyze
python3 scripts/controls_abc_workflow.py regression
