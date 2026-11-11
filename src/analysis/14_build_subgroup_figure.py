#!/usr/bin/env python3
"""Figure S12: all 54 subgroup effects on the same cohort-SD scale.

The 3 x 3 grid compares three metabolic indices across three mortality
outcomes. Points and 95% confidence intervals describe each age, sex and
diabetes stratum; the matching table reports full-cohort interaction tests.
"""
from pathlib import Path
import os, importlib.util
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

R = Path(os.environ['IR_REVISION_ROOT'])
C = R / '.'
F = Path(os.environ['IR_CR_STAGE']) / 'figures'
S = F / 'source_data'
S.mkdir(parents=True, exist_ok=True)
spec = importlib.util.spec_from_file_location('figure_style', R / 'src/figures/00_plot_helpers.py')
style = importlib.util.module_from_spec(spec)
spec.loader.exec_module(style)
plt.rcParams.update({'font.family': 'sans-serif', 'font.sans-serif': ['Arial', 'Helvetica', 'Liberation Sans'],
                     'pdf.fonttype': 42, 'svg.fonttype': 'none', 'savefig.dpi': 300,
                     'font.size': 7, 'axes.titlesize': 7, 'axes.labelsize': 7,
                     'xtick.labelsize': 6, 'ytick.labelsize': 6, 'figure.titlesize': 9})
sg = pd.read_csv(C / 'results/analysis/05_CR_M4_AO_M5_exploratory_subgroups.csv')
assert len(sg) == 54 and sg.exposure_scale.eq('per_sd').all()
assert (sg[['estimate', 'conf_low', 'conf_high']] > 0).all().all()
assert sg.groupby('exposure').index_sd.nunique().eq(1).all()
sg.to_csv(S / 'FigureS12_Subgroup_source.csv', index=False)
interactions = pd.read_csv(C / 'results/analysis/05_CR_M4_AO_subgroup_interactions.csv')
interactions.to_csv(S / 'FigureS10_Interaction_source.csv', index=False)
exposures = ['SHR', 'TyG', 'AIP']
outcomes = ['hospital_mortality', 'mortality_28d', 'icu_mortality']
titles = ['In-hospital death', '28-day death', 'First ICU death']
order = [('age', '<65'), ('age', '>=65'), ('sex', 'Female'), ('sex', 'Male'), ('diabetes', 'No'), ('diabetes', 'Yes')]
fig, axes = plt.subplots(3, 3, figsize=(183 / 25.4, 215 / 25.4))
for i, exposure in enumerate(exposures):
    for j, outcome in enumerate(outcomes):
        ax = axes[i, j]
        a = sg[(sg.exposure == exposure) & (sg.outcome == outcome)].set_index(['subgroup', 'level']).loc[order]
        positions = np.arange(6)[::-1]
        labels = [f'{dict(age="Age", sex="Sex", diabetes="Diabetes")[group]}: {level.replace(">=", "≥")}' for group, level in order]
        color = style.EXPOSURE_COLORS[exposure]
        ax.errorbar(a.estimate, positions, xerr=np.vstack([a.estimate-a.conf_low, a.conf_high-a.estimate]),
                    fmt='o', color=color, ecolor=color, ms=3, capsize=2, lw=1)
        ax.set_yticks(positions, labels, fontsize=5.5)
        ax.set_title(f'{chr(97+i*3+j)}  {exposure} · {titles[j]}', loc='left', fontsize=6.5, fontweight='bold')
        effect = 'HR' if outcome == 'mortality_28d' else 'sHR'
        style.style_forest_axis(ax, f'{effect} per cohort SD (95% CI)')
        style.set_readable_log_ticks(ax, min(1.0, a.conf_low.min()), max(1.0, a.conf_high.max()))
fig.suptitle('Subgroup associations · Model 3', fontweight='bold', y=.995)
fig.tight_layout(rect=(0, 0, 1, .97), w_pad=1.1, h_pad=1.8)
fig.canvas.draw()
# The author environment supplies the geometry auditor; exported data/figures
# remain reproducible without this optional external QA utility.
audit_root = os.environ.get('IR_FIGURE_QA_TOOLS')
if audit_root:
    import sys
    sys.path.insert(0, audit_root)
    from audit_panel_alignment import require_matplotlib_panel_alignment
    require_matplotlib_panel_alignment(fig, json_out=C / 'reports/internal/06_CR_subgroup_alignment.json')
fig.savefig(F / 'FigureS10_Subgroup_forest.pdf', bbox_inches='tight', dpi=300)
fig.savefig(F / 'FigureS10_Subgroup_forest.svg', bbox_inches='tight', dpi=300)
fig.savefig(F / 'FigureS10_Subgroup_forest.png', bbox_inches='tight', dpi=600)
plt.close(fig)
print('Updated Figure S10 from all 54 common-SD estimates and exported 27 interaction tests.')
