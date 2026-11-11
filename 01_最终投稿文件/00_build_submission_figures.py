#!/usr/bin/env python3
# Academic Figure Skill Asset Confirmation (verified against assets/figures/)
# (a) Cohort flow diagram -> cross-type inherit from SankeyDiagram -> param inherit
# (b) Forest plots -> cross-type inherit from BarComparison -> param inherit
# (c) ROC curves -> assets/figures/AUROC/plot_AUROC1.py -> visual adapt
# (d) RCS/KM curves -> assets/figures/LineTrend/plot_trend.py -> visual adapt
# (e) Incremental AUC -> assets/figures/BarComparison/plot_comparison_Trajectory.py -> visual adapt
# RULE: no production asset accepts the study-specific multi-panel data structure unchanged.

# Academic Figure Skill Typography Baseline — COPY VERBATIM, place at TOP of script
import matplotlib as mpl
mpl.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "Liberation Sans"],
    "font.size": 8,
    "axes.titlesize": 8,
    "axes.labelsize": 8,
    "xtick.labelsize": 7,
    "ytick.labelsize": 7,
    "legend.fontsize": 8,
    "figure.titlesize": 9,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "axes.linewidth": 0.6,
    "xtick.direction": "out",
    "ytick.direction": "out",
    "xtick.major.width": 0.6,
    "ytick.major.width": 0.6,
    "legend.frameon": False,
})

# Academic Figure Skill Nature/Cell/Science Color Palette -- COPY VERBATIM
CATEGORICAL = ["#2166AC", "#B2182B", "#1B7837", "#F1A340", "#762A83", "#666666"]
CATEGORICAL_EXTENDED = [
    "#2166AC", "#B2182B", "#1B7837", "#F1A340", "#762A83", "#666666",
    "#4393C3", "#D6604D", "#5AAE61", "#B35806", "#9970AB", "#999999",
]
DIVERGING   = ["#2166AC", "#F7F7F7", "#B2182B"]
SEQUENTIAL  = ["#F7FBFF", "#6BAED6", "#08306B"]
ACCENT_RED  = "#B2182B"
GREY        = "#999999"
BLACK       = "#222222"

# Academic Figure Skill Export Baseline — COPY VERBATIM
mpl.rcParams.update({
    "pdf.fonttype": 42,
    "svg.fonttype": "none",
    "savefig.bbox": "tight",
    "savefig.dpi": 300,
})

def save_cns_figure(fig, filename):
    """Standard Academic Figure Skill export: vector PDF + 300dpi PNG preview."""
    fig.savefig(f"{filename}.pdf", bbox_inches="tight", dpi=300)
    fig.savefig(f"{filename}.png", bbox_inches="tight", dpi=300)


from pathlib import Path
import math
import os
import shutil
import subprocess
import textwrap

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
from matplotlib.lines import Line2D
from matplotlib.patches import FancyBboxPatch
from docx import Document


REV = Path(os.environ.get("IR_REVISION_ROOT", Path(__file__).resolve().parents[1])).resolve()
ROOT = Path(os.environ.get("IR_PROJECT_ROOT", REV.parent)).resolve()
OUT = REV / "01_最终投稿文件"
FIG = OUT / "figures"
SRC = FIG / "source_data"
FIG.mkdir(parents=True, exist_ok=True)
SRC.mkdir(parents=True, exist_ok=True)

MODEL_COLORS = {"Model 1": "#4E79A7", "Model 2": "#F28E2B", "Model 3": "#2E8B57"}
MODEL_MARKERS = {"Model 1": "o", "Model 2": "s", "Model 3": "^"}
EXPOSURE_COLORS = {"SHR": "#2166AC", "TyG": "#F1A340", "AIP": "#762A83"}
TERTILE_COLORS = {"T1": "#2166AC", "T2": "#F1A340", "T3": "#B2182B"}
OUTCOME_ORDER = ["In-hospital mortality", "28-day mortality", "ICU mortality"]
EXPOSURE_ORDER = ["SHR", "TyG", "AIP"]


PATHS = {
    "M4_M5": REV / "02_PA_主分析_AdmissionOrigin_M5P/01_MIMICIV/03_汇总结果/01_PA_M4_AO_M5P_ALL_cox_results.csv",
    "M3_M5": REV / "02_PA_主分析_AdmissionOrigin_M5P/02_MIMICIII外部验证/03_汇总结果/02_PA_M3_AO_M5P_ALL_cox_results.csv",
    "NW_M5": REV / "02_PA_主分析_AdmissionOrigin_M5P/03_NWICU外部验证/03_汇总结果/03_PA_NW_AO_M5P_ALL_cox_results.csv",
    "M4_M20": REV / "04_SA_敏感性分析/01_AdmissionOrigin_M20P/01_MIMICIV/03_汇总结果/01_SA_M4_AO_M20P_ALL_cox_results.csv",
    "M3_M20": REV / "04_SA_敏感性分析/01_AdmissionOrigin_M20P/02_MIMICIII外部验证/03_汇总结果/02_SA_M3_AO_M20P_ALL_cox_results.csv",
    "NW_M20": REV / "04_SA_敏感性分析/01_AdmissionOrigin_M20P/03_NWICU外部验证/03_汇总结果/03_SA_NW_AO_M20P_ALL_cox_results.csv",
    "M4_L5": REV / "04_SA_敏感性分析/02_Landmark24h_M5P/01_MIMICIV/03_汇总结果/01_SA_M4_L24_M5P_ALL_cox_results.csv",
    "M3_L5": REV / "04_SA_敏感性分析/02_Landmark24h_M5P/02_MIMICIII外部验证/03_汇总结果/02_SA_M3_L24_M5P_ALL_cox_results.csv",
    "NW_L5": REV / "04_SA_敏感性分析/02_Landmark24h_M5P/03_NWICU外部验证/03_汇总结果/03_SA_NW_L24_M5P_ALL_cox_results.csv",
    "M4_L20": REV / "04_SA_敏感性分析/03_Landmark24h_M20P/01_MIMICIV/03_汇总结果/01_SA_M4_L24_M20P_ALL_cox_results.csv",
    "M3_L20": REV / "04_SA_敏感性分析/03_Landmark24h_M20P/02_MIMICIII外部验证/03_汇总结果/02_SA_M3_L24_M20P_ALL_cox_results.csv",
    "NW_L20": REV / "04_SA_敏感性分析/03_Landmark24h_M20P/03_NWICU外部验证/03_汇总结果/03_SA_NW_L24_M20P_ALL_cox_results.csv",
    "RCS": REV / "03_PSA_主分析补充分析/01_MIMICIV/03_汇总结果/01_PSA_M4_AO_M5P_ALL_rcs_source_data.csv",
    "KM": REV / "03_PSA_主分析补充分析/01_MIMICIV/03_汇总结果/01_PSA_M4_AO_NA_ALL_km_source_data.csv",
    "ROC_AUC": REV / "03_PSA_主分析补充分析/01_MIMICIV/03_汇总结果/01_PSA_M4_AO_NA_ALL_roc_auc.csv",
    "DELONG": REV / "03_PSA_主分析补充分析/01_MIMICIV/03_汇总结果/01_PSA_M4_AO_NA_ALL_delong.csv",
    "SCORE": REV / "03_PSA_主分析补充分析/01_MIMICIV/03_汇总结果/09_PSA_M4_AO_NA_ALL_score_comparison_cv_auc.csv",
    "COMPONENT": REV / "03_PSA_主分析补充分析/01_MIMICIV/03_汇总结果/08_PSA_M4_AO_NA_SHR_component_cv_auc.csv",
    "CALIBRATION": REV / "03_PSA_主分析补充分析/01_MIMICIV/03_汇总结果/10_PSA_M4_AO_M5P_ALL_calibration_brier_auc.csv",
    "DCA": REV / "03_PSA_主分析补充分析/01_MIMICIV/03_汇总结果/12_PSA_M4_AO_M5P_ALL_dca.csv",
    "SUBGROUP": REV / "03_PSA_主分析补充分析/01_MIMICIV/03_汇总结果/13_PSA_M4_AO_M5P_ALL_exploratory_subgroups.csv",
    "BASIC": REV / "03_PSA_主分析补充分析/01_MIMICIV/03_汇总结果/14_PSA_M4_AO_NA_ALL_basic_model_cv_auc.csv",
    "BASIC_NRI_IDI": REV / "03_PSA_主分析补充分析/01_MIMICIV/03_汇总结果/16_PSA_M4_AO_NA_ALL_basic_model_nri_idi.csv",
    "PER_SD": REV / "03_PSA_主分析补充分析/01_MIMICIV/03_汇总结果/01_PSA_M4_AO_M5P_ALL_per_sd_cox_results.csv",
    "LAB": REV / "04_SA_敏感性分析/05_实验室取值时间/01_MIMICIV/03_汇总结果/01_SA_M4_AO_M5P_ALL_lab_timing_results.csv",
    "LOGISTIC": REV / "04_SA_敏感性分析/06_FineGray与固定结局Logistic/01_MIMICIV/03_汇总结果/01_SA_M4_AO_M5P_ALL_logistic_results.csv",
    "FINEGRAY": REV / "04_SA_敏感性分析/06_FineGray与固定结局Logistic/01_MIMICIV/03_汇总结果/01_SA_M4_AO_M5P_ALL_finegray_results.csv",
    "WEIGHT": REV / "04_SA_敏感性分析/07_选择加权/01_MIMICIV/03_汇总结果/01_SA_M4_AO_M5P_ALL_selection_weighted_results.csv",
    "EXCLUSION": REV / "04_SA_敏感性分析/08_排除特定合并症/01_MIMICIV/03_汇总结果/01_SA_M4_AO_M5P_ALL_exclusion_results.csv",
}


def read(path):
    return pd.read_csv(path)


COHORT_UNITS = [
    ("01", "M4", "MIMIC-IV", "01_MIMICIV"),
    ("02", "M3", "MIMIC-III", "02_MIMICIII外部验证"),
    ("03", "NW", "NWICU", "03_NWICU外部验证"),
]


def collect_advanced(section, suffix):
    blocks = []
    for order, code, label, unit in COHORT_UNITS:
        path = REV / "04_SA_敏感性分析" / section / unit / "03_汇总结果" / f"{order}_SA_{code}_AO_M5P_ALL_{suffix}.csv"
        if path.exists():
            d = read(path); d["cohort_label"] = label; blocks.append(d)
    if not blocks:
        raise FileNotFoundError(f"No advanced sensitivity files for {section}/{suffix}")
    return pd.concat(blocks, ignore_index=True)


def normal_model(v):
    s = str(v).replace("Model", "Model ").replace("  ", " ").strip()
    return s


def normal_scale(v):
    return str(v).replace("-", "_").lower().strip()


def normal_outcome(v):
    return {
        "hospital_mortality": "In-hospital mortality",
        "hospital_outcome": "In-hospital mortality",
        "院内死亡": "In-hospital mortality",
        "mortality_28d": "28-day mortality",
        "outcome28": "28-day mortality",
        "28天死亡": "28-day mortality",
        "icu_mortality": "ICU mortality",
        "icu_outcome": "ICU mortality",
        "ICU死亡": "ICU mortality",
    }.get(str(v), str(v))


def clean_frame(d, scale="per_unit", exposures=EXPOSURE_ORDER):
    d = d.copy()
    d["model"] = d["model"].map(normal_model)
    d["outcome_label"] = d["outcome"].map(normal_outcome)
    d["scale_normalized"] = d["exposure_scale"].map(normal_scale)
    d = d[d["exposure"].isin(exposures) & (d["scale_normalized"] == scale)]
    for col in ["estimate", "conf_low", "conf_high", "p_value", "n", "events"]:
        d[col] = pd.to_numeric(d[col], errors="coerce")
    return d


def clean_results(path, scale="per_unit", exposures=EXPOSURE_ORDER):
    return clean_frame(read(path), scale=scale, exposures=exposures)


def p_text(p):
    if pd.isna(p):
        return "NA"
    return "<0.001" if p < 0.001 else f"{p:.3f}"


def validate_effects(d, name):
    required = ["estimate", "conf_low", "conf_high"]
    if d.empty:
        raise ValueError(f"{name}: no rows")
    if d[required].isna().any().any():
        raise ValueError(f"{name}: missing effect values")
    if not ((d["conf_low"] > 0) & (d["estimate"] > 0) & (d["conf_high"] > 0)).all():
        raise ValueError(f"{name}: non-positive estimate or CI")
    if not ((d["conf_low"] <= d["estimate"]) & (d["estimate"] <= d["conf_high"])).all():
        raise ValueError(f"{name}: estimate outside CI")


def model_legend():
    return [
        Line2D([0], [0], color=MODEL_COLORS[m], marker=MODEL_MARKERS[m], lw=1.3, ms=4, label=m)
        for m in ["Model 1", "Model 2", "Model 3"]
    ]


def style_forest_axis(ax, xlabel):
    ax.axvline(1, color="#777777", lw=0.7, ls="--", zorder=0)
    ax.set_xscale("log")
    ax.set_xlabel(xlabel)
    ax.grid(axis="x", color="#DDDDDD", lw=0.35, alpha=0.6)
    ax.spines["left"].set_visible(False)
    ax.tick_params(axis="y", length=0)


def set_readable_log_ticks(ax, low, high, max_ticks=3, min_ratio=1.35):
    """Use sparse, well-separated ticks for narrow publication forest panels."""
    low = max(float(low) * 0.88, 0.05)
    high = float(high) * 1.12
    candidates = np.array([
        0.10, 0.15, 0.20, 0.25, 0.30, 0.40, 0.50, 0.60, 0.75,
        1.00, 1.25, 1.50, 2.00, 2.50, 3.00, 4.00, 5.00, 6.00,
        8.00, 10.00, 15.00, 20.00,
    ])
    available = candidates[(candidates >= low) & (candidates <= high)]
    if len(available) >= 2:
        # Anchor the scale at its outer readable values. Add one central tick
        # only when all adjacent labels retain enough log-scale separation.
        ticks = [float(available[0]), float(available[-1])]
        if max_ticks >= 3 and ticks[-1] / ticks[0] >= min_ratio ** 2:
            target = math.sqrt(ticks[0] * ticks[-1])
            middle = sorted(
                (float(v) for v in available[1:-1]),
                key=lambda v: abs(math.log(v / target)),
            )
            for value in middle:
                if value / ticks[0] >= min_ratio and ticks[-1] / value >= min_ratio:
                    ticks.insert(1, value)
                    break
        ticks = np.asarray(ticks)
    elif len(available) == 1:
        ticks = available
    else:
        ticks = np.geomspace(low, high, 2)
    ax.set_xlim(low, high)
    ax.set_xticks(ticks)
    ax.xaxis.set_major_formatter(mticker.FuncFormatter(
        lambda value, _: (f"{value:.2f}".rstrip("0").rstrip("."))
    ))
    ax.xaxis.set_minor_locator(mticker.NullLocator())


def forest_by_exposure(d, basename, title, label_cols=("outcome_label",), effect_label="HR (95% CI)",
                       model_colors=True, model_col="model", category_order=None, legend_title=None,
                       width=183, height=145):
    validate_effects(d, basename)
    fig, axes = plt.subplots(1, 3, figsize=(width / 25.4, height / 25.4), sharex=False)
    offsets = {"Model 1": 0.22, "Model 2": 0.0, "Model 3": -0.22}
    fallback_vals = list(dict.fromkeys(d[model_col].astype(str)))
    fallback_colors = dict(zip(fallback_vals, CATEGORICAL_EXTENDED))
    for idx, exposure in enumerate(EXPOSURE_ORDER):
        ax = axes[idx]
        x = d[d["exposure"] == exposure].copy()
        x["category"] = x[list(label_cols)].astype(str).agg(" · ".join, axis=1)
        cats = category_order or list(dict.fromkeys(x["category"]))
        pos = {c: len(cats) - 1 - i for i, c in enumerate(cats)}
        for _, r in x.iterrows():
            key = str(r[model_col])
            off = offsets.get(key, (fallback_vals.index(key) - (len(fallback_vals) - 1) / 2) * 0.15)
            y = pos[r["category"]] + off
            color = MODEL_COLORS.get(key, fallback_colors[key]) if model_colors else fallback_colors[key]
            marker = MODEL_MARKERS.get(key, "o")
            ax.errorbar(
                r["estimate"], y,
                xerr=[[r["estimate"] - r["conf_low"]], [r["conf_high"] - r["estimate"]]],
                fmt=marker, color=color, ecolor=color, elinewidth=1.0, capsize=2.0, ms=3.5, zorder=3,
            )
        ax.set_yticks([pos[c] for c in cats])
        # All three exposure panels use the same row order. Printing the row
        # labels once preserves alignment and gives the effect axes enough
        # horizontal room in dense sensitivity figures.
        ax.set_yticklabels(cats if idx == 0 else [])
        ax.set_title(exposure, fontweight="bold")
        style_forest_axis(ax, effect_label)
        set_readable_log_ticks(ax, x["conf_low"].min(), x["conf_high"].max())
    handles = model_legend() if model_colors else [
        Line2D([0], [0], color=fallback_colors[k], marker="o", lw=1.2, label=k) for k in fallback_vals
    ]
    # Reserve separate top bands for the title and legend. This also prevents
    # tight-layout from expanding the canvas around figure-level artists.
    fig.tight_layout(rect=(0.0, 0.0, 1.0, 0.84), w_pad=1.15)
    fig.legend(handles=handles, loc="upper center", bbox_to_anchor=(0.5, 0.925), ncol=min(len(handles), 5),
               title=legend_title)
    fig.suptitle(title, y=0.985, fontweight="bold")
    save_cns_figure(fig, FIG / basename)
    plt.close(fig)


def build_flowchart():
    flow_source = SRC / "Figure1_detailed_cohort_flow_source.csv"
    flow = read(flow_source)
    required = {"cohort", "step_order", "criterion", "before_n", "excluded_n", "after_n", "stage"}
    missing = required.difference(flow.columns)
    if missing:
        raise ValueError(f"Figure 1 source is missing columns: {sorted(missing)}")
    for cohort, x in flow.groupby("cohort", sort=False):
        x = x.sort_values("step_order")
        if int(x.iloc[0]["before_n"]) != int(x.iloc[0]["after_n"]):
            raise ValueError(f"{cohort}: source row must have before_n = after_n")
        filtered = x[x["excluded_n"].notna()]
        if not (filtered["before_n"].astype(int) - filtered["excluded_n"].astype(int)
                == filtered["after_n"].astype(int)).all():
            raise ValueError(f"{cohort}: exclusion arithmetic failed")
        previous = None
        for _, row in x.iterrows():
            if previous is not None and int(row["before_n"]) != previous:
                raise ValueError(f"{cohort}: non-contiguous flow before step {row['step_order']}")
            previous = int(row["after_n"])

    fig, ax = plt.subplots(figsize=(183 / 25.4, 238 / 25.4))
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")

    def box(x, y, w, h, text, fc="#F7FBFF", ec="#4E79A7", size=5.6, weight="normal"):
        patch = FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.007,rounding_size=0.010",
                               facecolor=fc, edgecolor=ec, linewidth=0.8)
        ax.add_patch(patch)
        ax.text(x + w / 2, y + h / 2, text, ha="center", va="center", fontsize=size,
                fontweight=weight, linespacing=1.15)

    arrow = dict(arrowstyle="-|>", lw=0.75, color="#777777", shrinkA=3, shrinkB=3)
    labels = {"MIMIC-IV": "MIMIC-IV v3.1", "MIMIC-III": "MIMIC-III", "NWICU": "NWICU"}
    xs = [0.022, 0.347, 0.672]
    panel_w = 0.306
    top = 0.915
    bottom = 0.095
    ax.text(0.5, 0.982, "Detailed cohort selection", ha="center", va="top",
            fontsize=9, fontweight="bold", color=BLACK)
    for panel_x, cohort in zip(xs, ["MIMIC-IV", "MIMIC-III", "NWICU"]):
        x = flow[flow["cohort"] == cohort].sort_values("step_order").copy()
        ax.text(panel_x + panel_w / 2, 0.947, labels[cohort], ha="center", va="center",
                fontsize=7.2, fontweight="bold", color="#2166AC")
        n_rows = len(x)
        box_h = min(0.069, (top - bottom - 0.013 * (n_rows - 1)) / n_rows)
        gap = (top - bottom - n_rows * box_h) / max(n_rows - 1, 1)
        y_positions = [top - box_h - i * (box_h + gap) for i in range(n_rows)]
        for i, ((_, row), y) in enumerate(zip(x.iterrows(), y_positions)):
            stage = str(row["stage"])
            if stage == "source":
                fc, ec, size, weight = "#E8F1F8", "#2166AC", 5.9, "bold"
                text_value = f"{row['criterion']}\nN = {int(row['after_n']):,}"
            elif stage == "laboratory":
                fc, ec, size, weight = "#FFF8EE", "#F1A340", 5.15, "normal"
                criterion = str(row["criterion"]).replace(
                    " available and within data-quality range", "\nData-quality range"
                )
                text_value = (f"{criterion}\nExcluded n = {int(row['excluded_n']):,}; "
                              f"remaining N = {int(row['after_n']):,}")
            elif stage == "final":
                fc, ec, size, weight = "#EAF5EC", "#1B7837", 6.0, "bold"
                text_value = f"{row['criterion']}\nN = {int(row['after_n']):,}"
            else:
                fc, ec, size, weight = "#F7FBFF", "#4393C3", 5.35, "normal"
                text_value = (f"{row['criterion']}\nExcluded n = {int(row['excluded_n']):,}; "
                              f"remaining N = {int(row['after_n']):,}")
            box(panel_x, y, panel_w, box_h, text_value, fc=fc, ec=ec, size=size, weight=weight)
            if i < n_rows - 1:
                next_y = y_positions[i + 1]
                ax.annotate("", (panel_x + panel_w / 2, next_y + box_h),
                            (panel_x + panel_w / 2, y), arrowprops=arrow)

    ax.text(
        0.5, 0.041,
        "Laboratory exclusions are sequential in the displayed order (glucose → HbA1c → triglycerides → HDL-C); "
        "tests were assessed from 6 h before to 24 h after ICU admission. Exposure-defining values were not imputed.",
        ha="center", va="center", fontsize=5.2, color="#444444", wrap=True,
    )
    save_cns_figure(fig, FIG / "Figure1_Cohort_flowchart")
    plt.close(fig)


def build_primary_and_sensitivity_forests():
    main = clean_results(PATHS["M4_M5"])
    main.to_csv(SRC / "Figure2_primary_forest_source.csv", index=False)
    forest_by_exposure(
        main, "Figure2_Primary_M5_forest",
        "Primary MIMIC-IV admission-origin associations",
        effect_label="Hazard ratio (95% CI)", width=183, height=105,
    )

    groups = [
        ("FigureS01_Primary_M5_all_cohorts", "Primary admission-origin M5 Rubin-pooled results",
         ["M4_M5", "M3_M5", "NW_M5"]),
        ("FigureS02_Landmark24h_M5_all_cohorts", "24-hour landmark M5 Rubin-pooled results",
         ["M4_L5", "M3_L5", "NW_L5"]),
        ("FigureS03_AdmissionOrigin_M20_all_cohorts", "Admission-origin M20 Rubin-pooled results",
         ["M4_M20", "M3_M20", "NW_M20"]),
        ("FigureS04_Landmark24h_M20_all_cohorts", "24-hour landmark M20 Rubin-pooled results",
         ["M4_L20", "M3_L20", "NW_L20"]),
    ]
    for basename, title, keys in groups:
        blocks = []
        for key in keys:
            x = clean_results(PATHS[key])
            x["cohort_label"] = {"M4": "MIMIC-IV", "M3": "MIMIC-III", "NW": "NWICU"}[key.split("_")[0]]
            blocks.append(x)
        d = pd.concat(blocks, ignore_index=True)
        d.to_csv(SRC / f"{basename}_source.csv", index=False)
        cats = [f"{c} · {o}" for c in ["MIMIC-IV", "MIMIC-III", "NWICU"] for o in OUTCOME_ORDER]
        forest_by_exposure(d, basename, title, label_cols=("cohort_label", "outcome_label"),
                           category_order=cats, effect_label="Hazard ratio (95% CI)", width=183, height=205)

    blocks = []
    for cohort in ["M4", "M3", "NW"]:
        for imp in ["M5", "M20"]:
            x = clean_results(PATHS[f"{cohort}_{imp}"])
            x = x[x["model"] == "Model 3"].copy()
            x["cohort_label"] = {"M4": "MIMIC-IV", "M3": "MIMIC-III", "NW": "NWICU"}[cohort]
            x["imputation_label"] = imp
            blocks.append(x)
    comp = pd.concat(blocks, ignore_index=True)
    comp.to_csv(SRC / "FigureS05_M5_M20_comparison_source.csv", index=False)
    forest_by_exposure(comp, "FigureS05_M5_M20_comparison",
                       "M5 versus M20 admission-origin Model 3 estimates",
                       label_cols=("cohort_label", "outcome_label"), model_colors=False,
                       model_col="imputation_label",
                       category_order=[f"{c} · {o}" for c in ["MIMIC-IV", "MIMIC-III", "NWICU"] for o in OUTCOME_ORDER],
                       effect_label="Hazard ratio (95% CI)", legend_title="Imputation", width=183, height=200)


def build_rcs():
    d = read(PATHS["RCS"])
    d["outcome_label"] = d["outcome"].map(normal_outcome)
    d.to_csv(SRC / "Figure3_RCS_molar_scale_source.csv", index=False)
    fig, axes = plt.subplots(3, 3, figsize=(183 / 25.4, 220 / 25.4))
    for i, exposure in enumerate(EXPOSURE_ORDER):
        for j, outcome in enumerate(OUTCOME_ORDER):
            ax = axes[i, j]
            x = d[(d["exposure"] == exposure) & (d["outcome_label"] == outcome)].sort_values("x")
            ax.fill_between(x["x"].to_numpy(), x["conf_low"].to_numpy(), x["conf_high"].to_numpy(),
                            color=EXPOSURE_COLORS[exposure], alpha=0.18, lw=0)
            ax.plot(x["x"], x["estimate"], color=EXPOSURE_COLORS[exposure], lw=1.4)
            ax.axhline(1, color="#777777", lw=0.6, ls="--")
            ax.axvline(float(x["reference"].iloc[0]), color="#AAAAAA", lw=0.5, ls=":")
            ax.set_yscale("log")
            ax.set_title(f"{exposure} · {outcome}", fontweight="bold")
            ax.set_xlabel(exposure)
            ax.set_ylabel("Hazard ratio")
            ax.text(0.03, 0.97,
                    f"P-overall {p_text(float(x['p_overall'].iloc[0]))}\n"
                    f"P-nonlinear {p_text(float(x['p_nonlinear'].iloc[0]))}",
                    transform=ax.transAxes, va="top", fontsize=6.2)
    fig.suptitle("Restricted cubic spline associations (Model 3)", y=1.01, fontweight="bold")
    fig.tight_layout()
    save_cns_figure(fig, FIG / "Figure3_RCS")
    plt.close(fig)


def risk_at_times(x, times):
    out = []
    for tertile in ["T1", "T2", "T3"]:
        z = x[x["tertile"] == tertile].sort_values("time")
        vals = []
        initial = int(z["n_risk"].max()) if not z.empty else 0
        for t in times:
            eligible = z[z["time"] <= t]
            vals.append(initial if eligible.empty else int(eligible.iloc[-1]["n_risk"]))
        out.append(vals)
    return out


def build_km():
    d = read(PATHS["KM"])
    d["outcome_label"] = d["outcome"].map(normal_outcome)
    logrank = read(SRC / "Figure4_KM_logrank_summary.csv")
    logrank["outcome_label"] = logrank["outcome"].map(normal_outcome)
    risk = read(SRC / "Figure4_KM_fixed_time_risk_table.csv")
    risk["outcome_label"] = risk["outcome"].map(normal_outcome)
    fig, axes = plt.subplots(3, 3, figsize=(183 / 25.4, 240 / 25.4))
    for i, exposure in enumerate(EXPOSURE_ORDER):
        for j, outcome in enumerate(OUTCOME_ORDER):
            ax = axes[i, j]
            x = d[(d["exposure"] == exposure) & (d["outcome_label"] == outcome)].copy()
            for tertile in ["T1", "T2", "T3"]:
                z = x[x["tertile"] == tertile].sort_values("time")
                ax.step(z["time"], z["survival"], where="post", color=TERTILE_COLORS[tertile], lw=1.1,
                        label=tertile)
            rt = risk[(risk["exposure"] == exposure) & (risk["outcome_label"] == outcome)].copy()
            times = sorted(rt["time"].unique())
            risk_lines = [f"At risk (days {', '.join(str(int(t)) for t in times)})"]
            for tertile in ["T1", "T2", "T3"]:
                vals = rt[rt["tertile"] == tertile].sort_values("time")["n_risk"].astype(int).tolist()
                risk_lines.append(f"{tertile}: " + "  ".join(str(v) for v in vals))
            risk_text = "\n".join(risk_lines)
            lr = logrank[(logrank["exposure"] == exposure) & (logrank["outcome_label"] == outcome)].iloc[0]
            ax.set_title(f"{exposure} · {outcome}", fontweight="bold")
            ax.set_xlabel("Follow-up, days")
            ax.set_ylabel("Survival probability")
            ax.set_xlim(left=0)
            ax.set_ylim(max(0.45, float(x["survival"].min()) - 0.03), 1.005)
            ax.text(0.02, 0.03, f"Log-rank P {p_text(lr['logrank_p'])}", transform=ax.transAxes, fontsize=6)
            ax.text(0.5, -0.29, risk_text, transform=ax.transAxes,
                    ha="center", va="top", fontsize=4.8, linespacing=1.15)
    fig.legend(handles=[Line2D([0], [0], color=TERTILE_COLORS[t], lw=1.5, label=t) for t in ["T1", "T2", "T3"]],
               loc="upper center", ncol=3, bbox_to_anchor=(0.5, 1.005))
    fig.suptitle("Kaplan–Meier survival by index tertile (descriptive, unadjusted)", y=1.025, fontweight="bold")
    fig.subplots_adjust(hspace=0.86, wspace=0.42, top=0.94, bottom=0.08)
    save_cns_figure(fig, FIG / "Figure4_KM")
    plt.close(fig)


def build_roc():
    coords = read(SRC / "Figure5_ROC_curve_coordinates.csv")
    auc = read(PATHS["ROC_AUC"])
    auc["outcome_label"] = auc["outcome"].map(normal_outcome)
    coords["outcome_label"] = coords["outcome"].map(normal_outcome)
    fig, axes = plt.subplots(1, 3, figsize=(183 / 25.4, 67 / 25.4))
    for ax, outcome in zip(axes, OUTCOME_ORDER):
        for exposure in EXPOSURE_ORDER:
            x = coords[(coords["outcome_label"] == outcome) & (coords["exposure"] == exposure)]
            a = auc[(auc["outcome_label"] == outcome) & (auc["exposure"] == exposure)].iloc[0]
            ax.plot(x["false_positive_rate"], x["true_positive_rate"], color=EXPOSURE_COLORS[exposure],
                    lw=1.3, label=f"{exposure} (AUC {a['auc']:.3f})")
        ax.plot([0, 1], [0, 1], color="#999999", ls="--", lw=0.7)
        ax.set_xlim(0, 1)
        ax.set_ylim(0, 1)
        ax.set_aspect("equal", adjustable="box")
        ax.set_xlabel("False-positive rate")
        ax.set_ylabel("True-positive rate")
        ax.set_title(outcome, fontweight="bold")
        ax.legend(loc="lower right", fontsize=6)
    fig.suptitle("Standalone discrimination of metabolic indices", y=1.04, fontweight="bold")
    fig.tight_layout()
    save_cns_figure(fig, FIG / "Figure5_ROC")
    plt.close(fig)


def build_incremental_performance():
    bases = ["SOFA", "SAPS II", "APS III", "OASIS"]
    score_inc = read(PATHS["SCORE"]).copy()
    score_inc["outcome"] = score_inc["outcome"].map(normal_outcome)
    score_inc.to_csv(SRC / "Figure6_SeverityScore_increment_source.csv", index=False)
    basic_inc = read(PATHS["BASIC"]).copy()
    basic_inc["outcome"] = basic_inc["outcome"].map(normal_outcome)
    basic_inc.to_csv(SRC / "Figure6_BasicModel_increment_source.csv", index=False)
    read(PATHS["BASIC_NRI_IDI"]).to_csv(SRC / "TableS_BasicModel_NRI_IDI_source.csv", index=False)

    fig, axes = plt.subplots(2, 3, figsize=(183 / 25.4, 138 / 25.4))
    for j, outcome in enumerate(OUTCOME_ORDER):
        ax = axes[0, j]
        x = score_inc[score_inc["outcome"] == outcome]
        for k, exposure in enumerate(EXPOSURE_ORDER):
            z = x[x["exposure"] == exposure]
            xpos = np.arange(len(bases)) + (k - 1) * 0.13
            ax.scatter(xpos, z["extended_auc"], marker=["o", "s", "^"][k],
                       color=EXPOSURE_COLORS[exposure], s=14, label=exposure, zorder=3)
            for xx, (_, row) in zip(xpos, z.iterrows()):
                ax.plot([xx, xx], [row["baseline_auc"], row["extended_auc"]],
                        color=EXPOSURE_COLORS[exposure], lw=0.9)
                ax.scatter(xx, row["baseline_auc"], facecolor="white",
                           edgecolor=EXPOSURE_COLORS[exposure], s=11, lw=0.7, zorder=3)
        ax.set_xticks(range(len(bases)))
        ax.set_xticklabels(bases, rotation=30, ha="right")
        ax.set_ylabel("Cross-validated AUC")
        ax.set_title(outcome, fontweight="bold")
        ax.grid(axis="y", color="#DDDDDD", lw=0.35)

        ax2 = axes[1, j]
        z = basic_inc[basic_inc["outcome"] == outcome]
        base_auc = float(z["baseline_auc"].iloc[0])
        labels = ["Basic"] + [f"Basic+{e}" for e in EXPOSURE_ORDER]
        values = [base_auc] + [
            float(z[z["exposure"] == exposure]["extended_auc"].iloc[0])
            for exposure in EXPOSURE_ORDER
        ]
        colors = ["#777777"] + [EXPOSURE_COLORS[e] for e in EXPOSURE_ORDER]
        ax2.bar(np.arange(4), values, color=colors, width=0.64)
        for k, value in enumerate(values):
            label = f"{value:.3f}" if k == 0 else f"{value:.3f}\nΔ{value-base_auc:+.3f}"
            ax2.text(k, value + 0.003, label, ha="center", va="bottom", fontsize=5.0)
        ax2.set_xticks(range(4))
        ax2.set_xticklabels(labels, rotation=25, ha="right", fontsize=5.6)
        ax2.set_ylabel("Cross-validated AUC")
        ax2.set_title("Exploratory Basic Model increment", fontweight="bold")
        ax2.grid(axis="y", color="#DDDDDD", lw=0.35)
    fig.legend(handles=[Line2D([0], [0], color=EXPOSURE_COLORS[e], marker=["o", "s", "^"][i],
                               lw=1, label=e) for i, e in enumerate(EXPOSURE_ORDER)],
               loc="upper center", ncol=3, bbox_to_anchor=(0.5, 1.01))
    fig.suptitle("Incremental discrimination beyond severity scores and the exploratory Basic Model\n"
                 "Open points: baseline AUC; filled points: extended AUC; connecting segment: ΔAUC",
                 y=1.055, fontweight="bold")
    fig.tight_layout(h_pad=1.0)
    save_cns_figure(fig, FIG / "Figure6_Incremental_performance")
    plt.close(fig)


def build_subgroups():
    d = read(PATHS["SUBGROUP"])
    d["outcome_label"] = d["outcome"].map(normal_outcome)
    labels = {"age": "Age", "sex": "Sex", "diabetes": "Diabetes"}
    d["category"] = d.apply(lambda r: f"{labels.get(str(r['subgroup']).lower(), r['subgroup'])}: {r['level']}", axis=1)
    d.to_csv(SRC / "Figure7_Subgroup_forest_source.csv", index=False)
    fig, axes = plt.subplots(3, 3, figsize=(183 / 25.4, 210 / 25.4), sharex=False)
    for i, exposure in enumerate(EXPOSURE_ORDER):
        for j, outcome in enumerate(OUTCOME_ORDER):
            ax = axes[i, j]
            x = d[(d["exposure"] == exposure) & (d["outcome_label"] == outcome)].copy()
            cats = list(dict.fromkeys(x["category"]))
            pos = {c: len(cats) - 1 - k for k, c in enumerate(cats)}
            for _, r in x.iterrows():
                ax.errorbar(r["estimate"], pos[r["category"]],
                            xerr=[[r["estimate"] - r["conf_low"]], [r["conf_high"] - r["estimate"]]],
                            fmt="o", color=EXPOSURE_COLORS[exposure], ecolor=EXPOSURE_COLORS[exposure],
                            elinewidth=1, capsize=2, ms=3.5)
            ax.set_yticks([pos[c] for c in cats])
            ax.set_yticklabels(cats)
            ax.set_title(f"{exposure} · {outcome}", fontweight="bold")
            style_forest_axis(ax, "HR per SD (95% CI)")
            set_readable_log_ticks(ax, x["conf_low"].min(), x["conf_high"].max())
    fig.suptitle("Exploratory subgroup associations (Model 3; no multiplicity adjustment)",
                 y=1.01, fontweight="bold")
    fig.tight_layout()
    save_cns_figure(fig, FIG / "Figure7_Subgroup_forest")
    for old in [FIG / "FigureS11_Subgroup_forest.pdf", FIG / "FigureS11_Subgroup_forest.png"]:
        if old.exists():
            old.unlink()
    plt.close(fig)


def build_other_sensitivity_forests():
    lab = clean_frame(collect_advanced("05_实验室取值时间", "lab_timing_results"))
    lab["lab_label"] = lab["lab_definition"].replace({
        "first_available": "First available",
        "closest_to_admission": "Closest to ICU admission",
        "closest_to_icu": "Closest to ICU admission",
        "early_window_mean": "Early-window mean",
    })
    lab.to_csv(SRC / "FigureS06_Lab_timing_source.csv", index=False)
    cats = [f"{cohort} · {labdef} · {out}" for cohort in ["MIMIC-IV", "MIMIC-III", "NWICU"]
            for labdef in ["First available", "Closest to ICU admission", "Early-window mean"] for out in OUTCOME_ORDER]
    forest_by_exposure(lab, "FigureS06_Alternative_lab_timing", "Alternative laboratory timing",
                       label_cols=("cohort_label", "lab_label", "outcome_label"), category_order=cats,
                       effect_label="Hazard ratio (95% CI)", width=183, height=247)

    logistic = clean_frame(collect_advanced("06_FineGray与固定结局Logistic", "logistic_results"))
    logistic.to_csv(SRC / "FigureS07_Logistic_source.csv", index=False)
    forest_by_exposure(logistic, "FigureS07_Fixed_outcome_logistic", "Fixed-outcome logistic models",
                       label_cols=("cohort_label", "outcome_label"),
                       category_order=[f"{c} · {o}" for c in ["MIMIC-IV", "MIMIC-III", "NWICU"] for o in OUTCOME_ORDER],
                       effect_label="Odds ratio (95% CI)", width=183, height=190)

    fg = clean_frame(collect_advanced("06_FineGray与固定结局Logistic", "finegray_results"))
    fg.to_csv(SRC / "FigureS08_FineGray_source.csv", index=False)
    fg_order = [o for o in OUTCOME_ORDER if o in fg["outcome_label"].unique()]
    forest_by_exposure(fg, "FigureS08_FineGray", "Fine–Gray competing-risk models",
                       label_cols=("cohort_label", "outcome_label"),
                       category_order=[f"{c} · {o}" for c in ["MIMIC-IV", "MIMIC-III", "NWICU"] for o in fg_order],
                       effect_label="Subdistribution HR (95% CI)", width=183, height=155)

    w_all = clean_frame(collect_advanced("07_选择加权", "selection_weighted_results"), scale="per_sd")
    w_all.to_csv(SRC / "TableS_Selection_weighting_all_cohorts.csv", index=False)
    w = w_all[w_all["trimming"].astype(str).eq("p01_p99")].copy()
    w["trimming_label"] = "IPW trimmed at 1st–99th percentiles"
    w.to_csv(SRC / "FigureS09_Selection_weighting_source.csv", index=False)
    cats = [f"{c} · {out}" for c in ["MIMIC-IV", "MIMIC-III", "NWICU"] for out in OUTCOME_ORDER]
    forest_by_exposure(w, "FigureS09_Selection_weighting", "Selection-weighted sensitivity (IPW; 1st–99th percentile trimmed)",
                       label_cols=("cohort_label", "outcome_label"), category_order=cats,
                       model_colors=False, model_col="trimming_label",
                       effect_label="Weighted HR per SD", width=183, height=155)

    ex = clean_frame(collect_advanced("08_排除特定合并症", "exclusion_results"))
    ex["exclusion_label"] = ex["exclusion"].astype(str).str.replace("exclude_", "Exclude ").str.replace("_", " ")
    ex.to_csv(SRC / "TableS_Comorbidity_exclusions_all_cohorts.csv", index=False)
    for old in [FIG / "FigureS10_Comorbidity_exclusions.pdf", FIG / "FigureS10_Comorbidity_exclusions.png",
                SRC / "FigureS10_Comorbidity_exclusions_source.csv"]:
        if old.exists():
            old.unlink()

    tert = read(PATHS["M4_M5"]).copy()
    tert["model"] = tert["model"].map(normal_model)
    tert["outcome_label"] = tert["outcome"].map(normal_outcome)
    tert = tert[
        tert["exposure"].isin(EXPOSURE_ORDER)
        & tert["exposure_scale"].map(normal_scale).eq("tertile")
    ].copy()
    tert["tertile_label"] = tert["term"].replace({
        "T2相对T1": "T2 vs T1",
        "T3相对T1": "T3 vs T1",
    })
    for col in ["estimate", "conf_low", "conf_high", "p_value", "n", "events"]:
        tert[col] = pd.to_numeric(tert[col], errors="coerce")
    tert.to_csv(SRC / "FigureS11_Tertile_associations_source.csv", index=False)
    tert_cats = [
        f"{outcome} · {level}"
        for outcome in OUTCOME_ORDER
        for level in ["T2 vs T1", "T3 vs T1"]
    ]
    forest_by_exposure(
        tert, "FigureS11_Tertile_associations",
        "MIMIC-IV M5 Rubin-pooled tertile associations",
        label_cols=("outcome_label", "tertile_label"),
        category_order=tert_cats,
        effect_label="Hazard ratio (95% CI)",
        width=183, height=145,
    )

    per_sd = clean_results(PATHS["PER_SD"], scale="per_sd")
    per_sd.to_csv(SRC / "FigureS10_PerSD_associations_source.csv", index=False)
    forest_by_exposure(
        per_sd, "FigureS10_PerSD_associations",
        "MIMIC-IV M5 Rubin-pooled per-SD associations",
        effect_label="Hazard ratio per SD (95% CI)",
        width=183, height=105,
    )
    for old in [FIG / "FigureS12_PerSD_associations.pdf", FIG / "FigureS12_PerSD_associations.png",
                FIG / "FigureS13_Tertile_associations.pdf", FIG / "FigureS13_Tertile_associations.png",
                SRC / "FigureS12_PerSD_associations_source.csv", SRC / "FigureS13_Tertile_associations_source.csv"]:
        if old.exists():
            old.unlink()


def write_reports():
    stats = """# 图形统计与可复现性报告

- Figure 1：队列筛选计数来自锁定的纳入/排除结果；不使用患者级标识。
- Figure 2：MIMIC-IV admission-origin M5 Rubin-pooled Cox per-unit HR；每个点均由主分析 CSV 的 estimate、conf_low、conf_high 直接绘制。
- Figure 3：RCS 为 Model 3，5套插补模型在系数尺度按 Rubin 原则合并；曲线、95% CI、参考值、结点与 P-overall/P-nonlinear 均来自 RCS source CSV。
- Figure 4：KM 为未调整 tertile 描述性分析；曲线和风险人数来自 KM source CSV；log-rank P 从受限 MIMIC-IV 分析数据重新计算。
- Figure 5：ROC 坐标从受限 MIMIC-IV 共同样本重新计算；AUC 与现有锁定 AUC 汇总表核对，DeLong 比较来自现有 paired DeLong CSV。
- Figure 6：严重度评分与 Basic Model 增量均为重复10次5折交叉验证平均 AUC，并在每个基线—扩展模型对比内共用折分配；Basic Model 的变量筛选未嵌套于重采样，因此仅作探索性分析。其补充 NRI/IDI 为表观、post-selection 估计。
- Figure 7：探索性 subgroup Model 3 per-SD HR；未进行确认性多重比较校正。该图不在补充材料中重复。
- Supplementary Figure S9：对获得四项暴露定义实验室检查的概率建立稳定逆概率权重，并将权重截尾于第1和第99百分位；其目标是评估由已测变量解释的选择性检测偏倚。
- 合并症排除结果仅保留为补充数值表，不再制作信息密度过高的重复森林图。
- Supplementary Figures S1–S11：直接使用相应 M5、M20、landmark、laboratory timing、logistic、Fine–Gray、selection weighting、per-SD 和 tertile 汇总 CSV。

模型颜色固定：Model 1 蓝色 #4E79A7；Model 2 橙色 #F28E2B；Model 3 绿色 #2E8B57。所有森林图采用对数效应轴和无效线 1.0。
"""
    (FIG / "Figure_statistics_and_reproducibility_report.md").write_text(stats, encoding="utf-8")

    records = []
    for p in sorted(FIG.glob("*.png")):
        pdf = p.with_suffix(".pdf")
        records.append({"figure": p.stem, "png_exists": p.exists(), "pdf_exists": pdf.exists(),
                        "png_bytes": p.stat().st_size, "pdf_bytes": pdf.stat().st_size if pdf.exists() else 0})
    pd.DataFrame(records).to_csv(FIG / "figure_export_manifest.csv", index=False)


def main():
    subprocess.run(["Rscript", str(OUT / "00_export_roc_km_source.R")], check=True)
    build_flowchart()
    build_primary_and_sensitivity_forests()
    build_rcs()
    build_km()
    build_roc()
    build_incremental_performance()
    build_subgroups()
    build_other_sensitivity_forests()
    write_reports()
    print(f"Built {len(list(FIG.glob('*.png')))} PNG and {len(list(FIG.glob('*.pdf')))} PDF figure files")


if __name__ == "__main__":
    main()
