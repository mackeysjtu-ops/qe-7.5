#!/usr/bin/env python3
"""Generate the PNG figures for CAIMD_DEVELOPMENT_PROGRESS.html.

All time-series figures are drawn from raw evidence copies under
../assets/data (see recompute_metrics.py for parsing conventions).
No fictional data are drawn: every curve comes from a committed local
workflow file. Status/timeline figures visualize documented scalar facts.

Run with a Python that has numpy + matplotlib (no seaborn used), e.g.:
  PYTHONNOUSERSITE=1 ~/miniforge3/envs/pytorch/bin/python make_figures.py
"""
import json
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.dates as mdates  # noqa: F401  (kept for reference)
from matplotlib.patches import FancyBboxPatch
import numpy as np

from recompute_metrics import parse_eal, parse_caimd, parse_cel, det3, cond2, DT_PS

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "..", "assets", "data")
FIGS = os.path.join(HERE, "..", "assets", "figures")

# ----- palette (validated reference palette, light mode) --------------------
SURFACE = "#fcfcfb"
INK = "#0b0b0b"
INK2 = "#52514e"
MUTED = "#898781"
GRID = "#e1e0d9"
BASE = "#c3c2b7"
SERIES = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300"]
GOOD = "#0ca30c"
WARN = "#fab219"
CRIT = "#d03b3b"
BLUE = "#2a78d6"

plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Hiragino Sans GB", "Arial Unicode MS", "DejaVu Sans"],
    "font.size": 10.5,
    "figure.facecolor": SURFACE,
    "axes.facecolor": SURFACE,
    "savefig.facecolor": SURFACE,
    "axes.edgecolor": BASE,
    "axes.labelcolor": INK2,
    "xtick.color": MUTED,
    "ytick.color": MUTED,
    "axes.grid": True,
    "grid.color": GRID,
    "grid.linewidth": 0.8,
    "axes.axisbelow": True,
    "axes.unicode_minus": False,
    "figure.dpi": 100,
    "savefig.dpi": 200,
})


def style_ax(ax):
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    ax.grid(axis="x", visible=False)


def ascii_log_ticks(ax):
    """Plain-ASCII exponent tick labels (the CJK font lacks U+2212)."""
    from matplotlib.ticker import FuncFormatter, LogLocator
    ax.yaxis.set_major_locator(LogLocator(base=10.0))
    ax.yaxis.set_major_formatter(
        FuncFormatter(lambda v, p: "1e%d" % int(round(np.log10(v))) if v > 0 else ""))


def save(fig, name):
    path = os.path.join(FIGS, name)
    fig.savefig(path, bbox_inches="tight")
    plt.close(fig)
    print("wrote", path)


# ---------------------------------------------------------------------------
# Figure 3: development timeline (from git log dates + workflow reports)
# ---------------------------------------------------------------------------
def fig03_timeline():
    # dates verified against `git log --date=iso-strict` (assets/data/git_log.txt)
    # and the workflow directory timestamps / report generation times.
    events = [
        ("07-15", "a11913a", "QE 7.5 原始基线\nOriginal source baseline", "impl"),
        ("07-15", "f7256b1 → df2047e", "设计文档 + CAIMD V1 实现\nDesign doc + first strain-and-hold", "impl"),
        ("07-15", "7aa6191", "macOS 串行验证 16 例通过\nmacOS validation", "pass"),
        ("07-16", "d6dc6aa", "Leonardo GPU/MPI 验证 34 例通过\nA100 + MPI validation", "pass"),
        ("07-16", "77455a7 → 0282a6e", "边界加固 + 持久化验证 harness (29 例)\nHardening + tracked harness", "impl"),
        ("07-18", "b865b5d → 20a9f2f", "Ice X 100 GPa 示例框架\nIce X example framework", "impl"),
        ("07-18", "1f1f16d", "PW vc-md 显式速度输入支持\nExplicit ATOMIC_VELOCITIES path", "impl"),
        ("07-18", "13ba2eb", "300 K vc-md 初始化验证 (40 例)\nVelocity initialization validated", "pass"),
        ("07-22", "af7ebfc", "Nose thermostat + active CAIMD 支持\n(当前 HEAD)", "impl"),
        ("07-23", "本地 0.30 大应变", "投影精确但未掩码胞奇异\ncaimd projected cell is singular", "fail"),
        ("07-23", "本地 four-DOF", "Test1/Test2 升温至数千 K\nheating with/without free cell", "fail"),
        ("07-23", "本地 Control A", "普通 vc-md+Nose 无 CAIMD 也升温\nordinary baseline FAIL", "fail"),
        ("07-23", "Leonardo 50127909", "QE 7.4.1 拒绝 nose; af7ebfc 升温\ncomparison PARTIAL", "warn"),
    ]
    colors = {"impl": BLUE, "pass": GOOD, "fail": CRIT, "warn": WARN}
    labels = {"impl": "实现 / 加固 (commit)", "pass": "验证通过",
              "fail": "验证失败", "warn": "部分完成 / 受阻"}

    fig, ax = plt.subplots(figsize=(15.2, 6.2))
    n = len(events)
    ax.axhline(0, color=BASE, lw=1.2, zorder=1)
    heights = [1.0, -1.0, 2.05, -2.05]
    for i, (date, ref, desc, kind) in enumerate(events):
        h = heights[i % 4]
        ax.plot([i, i], [0, h], color=GRID, lw=1.0, zorder=1)
        ax.scatter([i], [0], s=64, color=colors[kind], zorder=3,
                   edgecolors=SURFACE, linewidths=1.5)
        va = "bottom" if h > 0 else "top"
        ax.text(i, h, "%s\n%s" % (ref, desc), ha="center", va=va,
                fontsize=7.6, color=INK,
                bbox=dict(boxstyle="round,pad=0.32", fc=SURFACE,
                          ec=colors[kind], lw=1.0))
        ax.text(i, 0.16 if h < 0 else -0.16, "2026-" + date, ha="center",
                va="bottom" if h < 0 else "top", fontsize=7.4, color=MUTED)
    handles = [plt.Line2D([], [], marker="o", ls="", color=colors[k],
                          label=labels[k]) for k in ("impl", "pass", "fail", "warn")]
    ax.legend(handles=handles, loc="lower center", frameon=False, fontsize=9,
              ncol=4, bbox_to_anchor=(0.5, -0.06))
    ax.set_xlim(-0.8, n - 0.2)
    ax.set_ylim(-3.5, 3.6)
    ax.axis("off")
    ax.set_title("CAIMD 开发与验证时间线（2026-07-15 至 2026-07-23，日期取自 git log 与工作流报告）",
                 color=INK, fontsize=12.5, pad=6)
    save(fig, "fig03_timeline.png")


# ---------------------------------------------------------------------------
# Figure 4: strain tracking in the failed high-rate 0.30 continuation run
# ---------------------------------------------------------------------------
def fig04_strain_tracking():
    c = parse_caimd(os.path.join(DATA, "cont_0p30_diamond8_continuous.caimd"))
    t = np.array(c["time_ps"])
    nom = np.array(c["nom_xx"])
    act = np.array(c["act_xx"])
    err = np.abs(act - nom)

    fig, (ax1, ax2) = plt.subplots(2, 1, sharex=True, figsize=(9.2, 6.4),
                                   gridspec_kw={"height_ratios": [2.0, 1.2]})
    ax1.plot(t, nom, color=INK2, lw=2.4, ls=(0, (5, 3)), label="nominal_xx（目标）")
    ax1.plot(t, act, color=BLUE, lw=1.6, label="actual_xx（实测）")
    ax1.set_ylabel("xx 应变（无量纲）")
    ax1.legend(frameon=False, loc="upper right")
    ax1.set_title("受控 xx 应变跟踪：本地 0.30 大应变连续加载（rate 0.05 ps$^{-1}$，失败前 265 行）",
                  color=INK, fontsize=11.5)
    style_ax(ax1)

    ax2.plot(t[err > 0], err[err > 0], color=BLUE, lw=1.4)
    ax2.set_yscale("log")
    ascii_log_ticks(ax2)
    ax2.set_ylabel("|actual - nominal|")
    ax2.set_xlabel("时间 time (ps)")
    ax2.annotate("最大偏差 max = %.2e（双精度舍入量级）" % err.max(),
                 xy=(t[np.argmax(err)], err.max()),
                 xytext=(0.30, 0.86), textcoords="axes fraction",
                 color=INK, fontsize=9.5,
                 arrowprops=dict(arrowstyle="-", color=MUTED, lw=0.9))
    style_ax(ax2)
    fig.align_ylabels()
    save(fig, "fig04_strain_tracking.png")


# ---------------------------------------------------------------------------
# Figure 5: temperature comparison across failed thermostat runs
# ---------------------------------------------------------------------------
def fig05_temperature_comparison():
    runs = [
        ("leonardo_af7ebfc.eal", "Leonardo af7ebfc 普通 vc-md+Nose（无 CAIMD）", SERIES[0]),
        ("control_A.eal", "本地 Control A 普通 vc-md+Nose（无 CAIMD）", SERIES[1]),
        ("test1.eal", "本地 Test1 active CAIMD 零载荷 + 自由 PR", SERIES[2]),
        ("test2.eal", "本地 Test2 active CAIMD 加载 + 全约束胞", SERIES[3]),
    ]
    fig, ax = plt.subplots(figsize=(11.6, 5.8))
    for fname, label, color in runs:
        rows = parse_eal(os.path.join(DATA, fname))
        t = np.array([s for s, _, _ in rows]) * DT_PS
        temp = np.array([x for _, _, x in rows])
        ax.plot(t, temp, color=color, lw=1.1, label=label, alpha=0.9)
    ax.axhline(300.0, color=MUTED, lw=1.4, ls=(0, (4, 3)))
    ax.text(3.1, 420, "300 K 目标 (tempw)", color=INK2, fontsize=9.2, ha="left")

    # first documented crossings of the Leonardo ordinary run (recomputed)
    leo = parse_eal(os.path.join(DATA, "leonardo_af7ebfc.eal"))
    lt = np.array([s for s, _, _ in leo]) * DT_PS
    ltemp = np.array([x for _, _, x in leo])
    for thr, txt, dy in ((600.0, "600 K @ 0.05 ps", 900), (2000.0, "2000 K @ 0.26 ps", 1400)):
        i = int(np.argmax(ltemp > thr))
        ax.scatter([lt[i]], [ltemp[i]], s=52, color=SERIES[0], zorder=5,
                   edgecolors=SURFACE, linewidths=1.2)
        ax.annotate(txt, xy=(lt[i], ltemp[i]), xytext=(lt[i] + 0.28, ltemp[i] + dy),
                    color=INK, fontsize=8.8,
                    arrowprops=dict(arrowstyle="-", color=MUTED, lw=0.9))
    t1 = parse_eal(os.path.join(DATA, "test1.eal"))
    ax.annotate("Test1", xy=(t1[-1][0] * DT_PS, t1[-1][2]),
                xytext=(t1[-1][0] * DT_PS + 0.06, t1[-1][2]),
                color=SERIES[2], fontsize=9, va="center")
    ax.set_xlabel("时间 time (ps)")
    ax.set_ylabel("离子动能温度 ionic temperature (K, 3N-3)")
    ax.set_title("温度对比：无论 CAIMD 是否激活，vc-md + Nose 轨迹均升温远超 300 K 目标（.eal 原始数据）",
                 color=INK, fontsize=11.5, pad=40)
    ax.legend(frameon=False, loc="lower center", fontsize=8.8, ncol=2,
              bbox_to_anchor=(0.5, 1.005))
    ax.set_xlim(0, 7.6)
    style_ax(ax)
    save(fig, "fig05_temperature_comparison.png")


# ---------------------------------------------------------------------------
# Figure 6: Leonardo af7ebfc — heating precedes cell deformation
# ---------------------------------------------------------------------------
def fig06_leonardo_temp_cell():
    leo = parse_eal(os.path.join(DATA, "leonardo_af7ebfc.eal"))
    t = np.array([s for s, _, _ in leo]) * DT_PS
    temp = np.array([x for _, _, x in leo])
    frames = parse_cel(os.path.join(DATA, "leonardo_af7ebfc.cel"))
    ct = np.array([s for s, _ in frames]) * DT_PS
    v0 = abs(det3(frames[0][1]))
    dv = np.array([100.0 * (abs(det3(m)) - v0) / v0 for _, m in frames])
    cond = np.array([cond2(m) for _, m in frames])

    fig, axes = plt.subplots(3, 1, sharex=True, figsize=(9.6, 7.6))
    axes[0].plot(t, temp, color=BLUE, lw=1.2)
    axes[0].axhline(300, color=MUTED, lw=1.2, ls=(0, (4, 3)))
    axes[0].set_ylabel("温度 T (K)")
    axes[0].set_title("Leonardo af7ebfc 普通 vc-md + Nose（无 CAIMD，1000 步）：升温先于任何显著胞变形",
                      color=INK, fontsize=11.5)
    axes[1].plot(ct, dv, color=BLUE, lw=1.2)
    axes[1].set_ylabel("体积变化 ΔV/V$_0$ (%)")
    axes[1].annotate("终态仅 +0.30%", xy=(ct[-1], dv[-1]), xytext=(1.44, 0.10),
                     color=INK, fontsize=9,
                     arrowprops=dict(arrowstyle="-", color=MUTED, lw=0.9))
    axes[2].plot(ct, cond, color=BLUE, lw=1.2)
    axes[2].set_ylabel("条件数 cond(H)")
    axes[2].set_xlabel("时间 time (ps)")
    axes[2].annotate("最大 1.017（远离奇异）", xy=(ct[int(np.argmax(cond))], cond.max()),
                     xytext=(0.62, 0.35), textcoords="axes fraction",
                     color=INK, fontsize=9,
                     arrowprops=dict(arrowstyle="-", color=MUTED, lw=0.9))
    for ax in axes:
        style_ax(ax)
        for x0, lab in ((0.05, "600 K"), (0.258, "2000 K")):
            ax.axvline(x0, color=WARN, lw=1.0, ls=(0, (2, 2)))
    axes[0].text(0.05, 8300, " 600 K", color=INK2, fontsize=8.4)
    axes[0].text(0.258, 8300, " 2000 K", color=INK2, fontsize=8.4)
    fig.align_ylabels()
    save(fig, "fig06_leonardo_temp_cell.png")


# ---------------------------------------------------------------------------
# Figure 7: validation status matrix
# ---------------------------------------------------------------------------
def fig07_status_matrix():
    rows = [
        ("输入接口 input interface (&CELL caimd*)", "PASS", "16/29/34/40 例 harness；解析/拒绝路径全覆盖"),
        ("应变事件调度 event scheduling (stride)", "PASS", "stride-2 事件序列 T,F,T,F,…；事件计数精确"),
        ("掩码应变投影 masked strain projection", "PASS", "nominal-actual ≤ 2.3e-16（0.30 失败轨迹中仍精确）"),
        ("未掩码 PR 响应 unmasked PR response", "PARTIAL", "短程正常；高速率大应变下演化至奇异胞"),
        ("反对称胞运动抑制 rotation suppression", "PASS", "‖W‖$_F$ ≤ 3.8e-08（Leonardo GPU 两进程）"),
        ("物理速度保持 Cartesian velocity preservation", "PASS", "变换残差 ≤ 1.6e-14（含失败轨迹）"),
        ("禁用/全假掩码 no-op regression", "PASS", "与基线逐位一致（能量/力/应力/轨迹）"),
        ("重启 restart", "NOT IMPLEMENTED", "V1 明确拒绝；需持久化 CAIMD+Nose 状态"),
        ("thermostat 输入路径 (nose allowlist)", "PASS", "not_controlled/nose 白名单；其余显式拒绝"),
        ("短程 active CAIMD + Nose 运行", "PASS", "本地定向 harness 通过（短 smoke）"),
        ("长时间温度稳定性 long-time Nose stability", "FAIL", "普通基线即升温至 ~8000 K（Control A / Leonardo）"),
        ("生产就绪 production readiness", "FAIL", "NOT READY：基线与延拓语义未解决"),
    ]
    colors = {"PASS": GOOD, "PARTIAL": WARN, "FAIL": CRIT,
              "NOT IMPLEMENTED": MUTED, "NOT YET TESTED": MUTED}
    fig, ax = plt.subplots(figsize=(12.4, 6.6))
    n = len(rows)
    for i, (label, status, note) in enumerate(rows):
        y = n - 1 - i
        ax.text(0.0, y, label, ha="left", va="center", fontsize=9.6, color=INK)
        c = colors[status]
        ax.add_patch(FancyBboxPatch((0.47, y - 0.30), 0.135, 0.60,
                                    boxstyle="round,pad=0.012,rounding_size=0.06",
                                    fc=c, ec="none", mutation_aspect=0.5))
        ax.text(0.5375, y, status, ha="center", va="center", fontsize=8.4,
                color="#ffffff", fontweight="bold")
        ax.text(0.625, y, note, ha="left", va="center", fontsize=8.6, color=INK2)
        if i < n - 1:
            ax.axhline(y - 0.5, color=GRID, lw=0.7)
    ax.set_xlim(-0.01, 1.18)
    ax.set_ylim(-0.6, n - 0.4)
    ax.axis("off")
    ax.set_title("CAIMD 验证状态矩阵（截至 2026-07-23，af7ebfc）", color=INK,
                 fontsize=12.5, pad=14, loc="left")
    save(fig, "fig07_status_matrix.png")


# ---------------------------------------------------------------------------
# Figure 9: all six actual strain components in the failed 0.30 run
# ---------------------------------------------------------------------------
def fig09_unmasked_components():
    c = parse_caimd(os.path.join(DATA, "cont_0p30_diamond8_continuous.caimd"))
    t = np.array(c["time_ps"])
    comps = ["xx", "yy", "zz", "yz", "xz", "xy"]
    fig, ax = plt.subplots(figsize=(10.4, 5.6))
    for i, k in enumerate(comps):
        y = np.array(c["act_" + k])
        lw = 2.4 if k == "xx" else 1.4
        ax.plot(t, y, color=SERIES[i], lw=lw, label=k)
    ax.axhline(0, color=BASE, lw=1.0)
    for k, dy in (("xz", 0.0), ("zz", 0.0), ("yy", 0.0)):
        i = comps.index(k)
        y = np.array(c["act_" + k])
        ax.annotate(k, xy=(t[-1], y[-1]), xytext=(t[-1] + 0.006, y[-1] + dy),
                    color=SERIES[i], fontsize=9.4, va="center")
    ax.annotate("掩码分量 xx 保持 -0.0265", xy=(t[-1], c["act_xx"][-1]),
                xytext=(0.30, 0.30), textcoords="axes fraction", color=SERIES[0],
                fontsize=9.6,
                arrowprops=dict(arrowstyle="-", color=MUTED, lw=0.9))
    ax.set_xlabel("时间 time (ps)")
    ax.set_ylabel("实际应变分量 actual strain (tensor)")
    ax.set_title("0.30 大应变失败轨迹：掩码 xx 分量保持精确，未掩码分量在 ~0.4 ps 后失稳发散",
                 color=INK, fontsize=11.5)
    ax.legend(frameon=False, loc="upper left", fontsize=9, ncol=6, title="分量 component",
              title_fontsize=8.5)
    style_ax(ax)
    save(fig, "fig09_unmasked_strain_components.png")


# ---------------------------------------------------------------------------
# Figure 10: velocity transformation residuals stay at roundoff scale
# ---------------------------------------------------------------------------
def fig10_velocity_jump():
    c = parse_caimd(os.path.join(DATA, "cont_0p30_diamond8_continuous.caimd"))
    t = np.array(c["time_ps"])
    vj = np.array(c["velocity_jump_max"])
    m = vj > 0
    fig, ax = plt.subplots(figsize=(9.6, 4.6))
    ax.plot(t[m], vj[m], color=BLUE, lw=1.2)
    ax.set_yscale("log")
    ascii_log_ticks(ax)
    ax.set_xlabel("时间 time (ps)")
    ax.set_ylabel("velocity_jump_max（a.u.，对数轴）")
    ax.set_title("速度基底变换残差：即使在失败的 0.30 轨迹中也保持在浮点舍入量级（≤ 1.6e-14）",
                 color=INK, fontsize=11.5)
    ax.annotate("max = %.2e" % vj.max(), xy=(t[np.argmax(vj)], vj.max()),
                xytext=(0.32, 0.88), textcoords="axes fraction", color=INK,
                fontsize=9.6,
                arrowprops=dict(arrowstyle="-", color=MUTED, lw=0.9))
    style_ax(ax)
    save(fig, "fig10_velocity_jump.png")


# ---------------------------------------------------------------------------
# Figure 11: Control A Nose extended-energy-like drift (from ABC analysis CSV)
# ---------------------------------------------------------------------------
def fig11_nose_extended():
    src = os.path.join(DATA, "controls_ABC_nose_extended_energy_timeseries.csv")
    import csv
    ts, tot = [], []
    with open(src) as fh:
        for row in csv.DictReader(fh):
            if row["control"] != "A":
                continue
            ts.append(float(row["time_ps"]))
            tot.append(float(row["extended_total_like_Ry"]))
    fig, ax = plt.subplots(figsize=(9.6, 4.6))
    ax.plot(ts, tot, color=BLUE, lw=1.3)
    ax.set_xlabel("时间 time (ps)")
    ax.set_ylabel("扩展总能类量 extended-total-like (Ry)")
    ax.set_title("Control A（普通 vc-md + Nose，无 CAIMD）：扩展能量类量随升温漂移，未见守恒平台",
                 color=INK, fontsize=11.5)
    style_ax(ax)
    save(fig, "fig11_nose_extended_energy.png")


if __name__ == "__main__":
    os.makedirs(FIGS, exist_ok=True)
    fig03_timeline()
    fig04_strain_tracking()
    fig05_temperature_comparison()
    fig06_leonardo_temp_cell()
    fig07_status_matrix()
    fig09_unmasked_components()
    fig10_velocity_jump()
    fig11_nose_extended()
    print("all figures written")
