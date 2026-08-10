# CAIMD Development Progress HTML Report

本目录包含 QE 7.5 CAIMD（constrained ab initio molecular dynamics）开发与验证
进度的可视化 HTML 报告及其全部本地证据副本。

## 主交付物

- `CAIMD_DEVELOPMENT_PROGRESS.html` — 中文进度报告（17 节，11 张图）。
  直接在浏览器中打开即可（双击或 `open CAIMD_DEVELOPMENT_PROGRESS.html`）。
  页面完全离线自包含：内联 CSS/JS，无外部字体、无 CDN、无网络依赖；
  图片通过相对路径引用 `assets/`，因此**整个目录**可作为一个单元移动或打包。

## 报告如何生成

1. 检查仓库状态（`main` @ `af7ebfc757c32e0b4f050ed9f2208e01020b2ee3`，
   基线 `a11913ac6cd37f39cf31056e8e4823ec7a5de49c`）与全部 CAIMD 文档。
2. 对照 af7ebfc 源码核实实现细节（`caimd_module.f90`、`vcsmd.f90`、
   `vcsubs.f90`、`input.f90`、`read_namelists.f90` 等）。
3. 将四个验证工作流的原始机器可读证据复制到 `assets/data/`
   （`.eal` / `.cel` / `.caimd` / CSV / JSON），关键报告复制到 `evidence/`。
4. `scripts/recompute_metrics.py` 从原始数据重算关键指标并与既有记录交叉核对
   （输出 `assets/data/recomputed_metrics.json`）。
5. `scripts/make_figures.py` 生成全部 matplotlib PNG 图；
   三张架构 / 概念 / 决策树图为手绘 SVG（`assets/architecture/`）。
6. 组装 HTML，执行链接 / 结构 / 自包含性校验，计算 SHA256
   （见 `REPORT_METADATA.json`）。

## 如何重新生成图

```bash
cd scripts
# 需要带 numpy + matplotlib (+ pandas 可选) 的 Python；本机可用：
PYTHONNOUSERSITE=1 ~/miniforge3/envs/pytorch/bin/python recompute_metrics.py
PYTHONNOUSERSITE=1 ~/miniforge3/envs/pytorch/bin/python make_figures.py
```

脚本只读取 `../assets/data/` 内的副本，不访问网络，不修改任何仓库文件。

## 底层证据位置

- 本目录副本：`assets/data/`（原始数据）、`evidence/`（报告与输入）。
- 原始工作流（未修改，绝对路径）：
  - `examples/caimd_diamond8_local_nose_continuation_0p30_20260723T102331Z`
  - `examples/caimd_diamond8_local_nose_four_dof_tests_0p03_20260723T115714Z`
  - `examples/caimd_diamond8_local_nose_controls_ABC_20260723T134303Z`
  - `examples/diamond8_ordinary_nose_qe74_vs_af7ebfc_completed_20260723T185307Z`
  - `examples/diamond8_ordinary_nose_qe74_vs_qe75_af7ebfc_20260723T175625Z`

## 保证

- 未修改任何已跟踪 QE 源文件（生成后经 `git status --short` 复核）。
- 未修改任何原始计算输出或既有报告；所有新文件仅在本目录内。
- 未使用 Leonardo / SSH / SCP / Slurm / 网络；未创建 Git 提交、未 push。
- 所有图均来自真实本地数据或文档化标量指标；无虚构曲线。
- 本地仓库中不存在 64 原子 Diamond 工作流证据；所有长时热浴测试均为
  8 原子 Diamond8 胞（报告 §8 有明确说明）。
