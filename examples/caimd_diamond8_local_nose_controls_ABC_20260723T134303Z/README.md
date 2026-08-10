# Diamond8 CAIMD Nose Controls A/B/C

Generated: 2026-07-23T13:53:34Z

Strictly local macOS diagnostic workflow for QE 7.5 CAIMD Nose behavior.

Controls:

- A: ordinary PW `vc-md`, Nose, CAIMD disabled.
- B: `caimd=.true.` all six masks false, zero increments, intended no-op.
- C: active all-true zero-strain CAIMD projection.

Run with:

```bash
python3 scripts/controls_abc_workflow.py prepare
python3 scripts/controls_abc_workflow.py run-all
python3 scripts/controls_abc_workflow.py analyze
python3 scripts/controls_abc_workflow.py regression
```

All generated files are inside this directory. The QE source tree is not
modified by this workflow.
