#!/bin/bash
set -euo pipefail
cd "$(dirname "$0")/../.."
python3 scripts/analyze_caimd_stress_strain.py
