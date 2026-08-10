#!/usr/bin/env python3
"""Aggregate all stage analysis_summary.json into FINAL_CHAIN_REPORT.md +
analysis/overall_chain_summary.json. Run at the end of JOB 4."""
import json
import os
import sys

ROOT = sys.argv[1] if len(sys.argv) > 1 else "."


def load(p):
    try:
        return json.load(open(p))
    except Exception:
        return None


def g(d, *keys, default=None):
    for k in keys:
        if not isinstance(d, dict):
            return default
        d = d.get(k, default)
    return d


vc = load(os.path.join(ROOT, "01_vc_relax", "analysis_summary.json"))
nvt = load(os.path.join(ROOT, "02_fixed_nvt", "analysis_summary.json"))
npt = load(os.path.join(ROOT, "03_free_npt", "analysis_summary.json"))
z = load(os.path.join(ROOT, "04_caimd", "zero_strain", "analysis_summary.json"))
nz = load(os.path.join(ROOT, "04_caimd", "small_nonzero_strain", "analysis_summary.json"))

# determine first failure / overall
overall = "FULL_CHAIN_VALIDATED"
first_fail = None
if not (vc and vc.get("gate_pass")):
    overall, first_fail = "VCRELAX_FAILED", "vc-relax"
elif not (nvt and nvt.get("gate_pass")):
    overall, first_fail = "NVT_TEMPERATURE_FAILED", "nvt"
elif not (npt and npt.get("gate_pass")):
    overall, first_fail = "NPT_TEMPERATURE_OR_PRESSURE_FAILED", "npt"
elif not (z and z.get("gate_pass")):
    overall, first_fail = "CAIMD_ZERO_STRAIN_FAILED", "caimd_zero"
elif not (nz and nz.get("gate_pass")):
    overall, first_fail = "CAIMD_NONZERO_STRAIN_FAILED", "caimd_nonzero"

summ = {"overall_status": overall, "first_failed_stage": first_fail,
        "stages": {"vc_relax": vc, "nvt": nvt, "npt": npt, "caimd_zero": z, "caimd_nonzero": nz},
        "no_qe_source_modified": True}
os.makedirs(os.path.join(ROOT, "analysis"), exist_ok=True)
json.dump(summ, open(os.path.join(ROOT, "analysis", "overall_chain_summary.json"), "w"), indent=2, default=float)


def f(x, n=1):
    try:
        return f"{float(x):.{n}f}"
    except (TypeError, ValueError):
        return "n/a"


L = []
L.append("# Final Chain Report: diamond8 vc-relax / NVT / NPT / CAIMD")
L.append("")
L.append(f"**Overall status: `{overall}`**" + (f" (first failure: {first_fail})" if first_fail else ""))
L.append("")
L.append("| # | stage | gate | key result |")
L.append("|---|---|---|---|")
L.append(f"| 1 | vc-relax | {g(vc,'gate_pass')} | final P {f(g(vc,'final_pressure_kbar'),2)} kbar, celldm(1) {f(g(vc,'celldm1_bohr'),4)} bohr |")
L.append(f"| 2 | fixed NVT | {g(nvt,'gate_pass')} | final-2ps T {f(g(nvt,'temperature','final_2ps_mean'))} K, P mean {f(g(nvt,'pressure','mean'),1)} kbar |")
L.append(f"| 3 | free NPT | {g(npt,'gate_pass')} | final-2ps T {f(g(npt,'temperature','final_2ps_mean'))} K, P block {f(g(npt,'pressure','final_block_mean_kbar'),1)} kbar, dV {f(g(npt,'cell','volume_change_vs_refvol_pct'),2)}% |")
L.append(f"| 4A | CAIMD zero | {g(z,'gate_pass')} | strain err {g(z,'caimd','max_abs_nominal_minus_actual')}, final-1.5ps T {f(g(z,'temperature','final_1p5ps_mean'))} K |")
L.append(f"| 4B | CAIMD nonzero | {g(nz,'gate_pass')} | events {g(nz,'caimd','stage4B_event_count')}, final nom xx {g(nz,'caimd','stage4B_final_nominal_xx')}, track err {g(nz,'caimd','stage4B_xx_tracking_error')} |")
L.append("")
L.append("## Required answers")
L.append(f"1. vc-relax reached 0 kbar: {g(vc,'gate_pass')} (final P {f(g(vc,'final_pressure_kbar'),2)} kbar, |P|<=5 kbar gate).")
L.append(f"2. Fixed-cell NVT maintained 300 K: {g(nvt,'gate_pass')} (final-2ps {f(g(nvt,'temperature','final_2ps_mean'))} K, block {f(g(nvt,'temperature','block_uncertainty','block_mean'))}+/-{f(g(nvt,'temperature','block_uncertainty','sem'))}).")
L.append(f"3. Fixed-cell NVT pressure: {f(g(nvt,'pressure','mean'),1)} kbar mean (diagnostic).")
L.append(f"4. Free-cell NPT maintained 300 K: {g(npt,'gate_pass') and (250 <= (g(npt,'temperature','final_2ps_mean') or 0) <= 350)} (final-2ps {f(g(npt,'temperature','final_2ps_mean'))} K).")
L.append(f"5. Free-cell NPT maintained 0 kbar: block mean {f(g(npt,'pressure','final_block_mean_kbar'),1)}+/-{f(g(npt,'pressure','block_sem_kbar'),1)} kbar.")
L.append(f"6. Ionic Nose active: NVT {g(nvt,'nose_ion_active')}, NPT {g(npt,'nose_ion_active')}.")
L.append(f"7. Cell Nose active (NPT): requested cell_temperature='nose' (temph=300, fnoseh=1); see NPT report.")
L.append(f"8. CAIMD disabled in JOB 3: {g(npt,'caimd_disabled')}.")
L.append(f"9. CAIMD zero strain stable: {g(z,'gate_pass')} (max strain err {g(z,'caimd','max_abs_nominal_minus_actual')}, vjump {g(z,'caimd','max_velocity_jump_max')}).")
L.append(f"10. Small nonzero strain applied correctly: {g(nz,'gate_pass')} (final nominal xx {g(nz,'caimd','stage4B_final_nominal_xx')}, expected cumulative {g(nz,'caimd','stage4B_expected_cumulative_xx')}, tracking err {g(nz,'caimd','stage4B_xx_tracking_error')}).")
L.append(f"11. First failed stage: {first_fail or 'none'}.")
L.append(f"12. Ready for larger/longer test: {'yes' if overall == 'FULL_CHAIN_VALIDATED' else 'no - resolve ' + str(first_fail) + ' first'}.")
L.append("13. QE source modified: No.")
L.append("")
L.append("Caveats: cutoffs 80/320 Ry (NC LDA C.UPF, no UPF recommendation); .tv/stdout T is a running "
         "average (instantaneous T uses 2*eka/((3N-3)k_B)); Nose chain energy not written -> extended "
         "energy is a partial diagnostic; 8-atom fluctuations are large - stability judged by block means.")
open(os.path.join(ROOT, "FINAL_CHAIN_REPORT.md"), "w").write("\n".join(L) + "\n")
print("overall:", overall)
