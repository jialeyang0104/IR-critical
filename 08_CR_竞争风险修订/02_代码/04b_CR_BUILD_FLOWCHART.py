#!/usr/bin/env python3
"""Draw the three-cohort selection flow with sequential side exclusions."""
from pathlib import Path
import os
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch

R = Path(os.environ['IR_REVISION_ROOT'])
C = R / '08_CR_竞争风险修订'
F = Path(os.environ['IR_CR_STAGE']) / 'figures'
S = F / 'source_data'
S.mkdir(parents=True, exist_ok=True)
flow = pd.read_csv(C / '03_汇总结果/00_CR_patient_flow.csv')
assert flow.cohort.tolist() == ['MIMIC-IV', 'MIMIC-III', 'NWICU']
flow.to_csv(S / 'Figure1_patient_flow_source.csv', index=False)

plt.rcParams.update({'font.family': 'sans-serif',
                     'font.sans-serif': ['Arial', 'Helvetica', 'Liberation Sans'],
                     'pdf.fonttype': 42, 'svg.fonttype': 'none',
                     'font.size': 6.7, 'savefig.dpi': 300})
fig = plt.figure(figsize=(183 / 25.4, 215 / 25.4))
ax = fig.add_axes([0, 0, 1, 1])
ax.set(xlim=(0, 183), ylim=(0, 215))
ax.axis('off')

def box(x, y, w, h, text, fill='white', bold=False, size=6.7):
    patch = FancyBboxPatch((x, y-h/2), w, h,
                           boxstyle='round,pad=0,rounding_size=1.2',
                           facecolor=fill, edgecolor='#202020', linewidth=.85)
    ax.add_patch(patch)
    label = ax.text(x+w/2, y, text, ha='center', va='center',
                    fontsize=size, color='#111111',
                    fontweight='bold' if bold else 'normal', linespacing=1.25)
    return patch, label

def arrow(a, b):
    ax.add_patch(FancyArrowPatch(a, b, arrowstyle='->',
                                 mutation_scale=8, linewidth=.85,
                                 color='#202020', shrinkA=0, shrinkB=0))

orange = '#F7CFAD'
blue = '#B8C4D3'
box(2, 207, 57, 10, 'Primary cohort', orange, True, 8)
box(63, 207, 118, 10, 'Association replication cohorts', blue, True, 8)
ys = [183, 151, 119, 87, 55, 23]
steps = []
text_boxes = []
for k, (_, row) in enumerate(flow.iterrows()):
    left = 2 + k*61
    fill = '#FBE8D8' if k == 0 else '#E2E7EE'
    ax.text(left+28.5, 195, row.cohort, ha='center', va='center',
            fontsize=8, fontweight='bold', color='#111111')
    labels = [
        f'All recorded ICU stays\nn = {row.all_icu_stays:,}',
        f'First ICU stay\nper patient\nn = {row.first_patient_stays:,}',
        f'Source age available\nn = {row.source_covariates_available:,}',
        f'Adults aged ≥18 years\nn = {row.adult_first_stays:,}',
        f'ICU stay and hospital\nfollow-up ≥24 h\nn = {row.eligible_pre_laboratory:,}',
        f'Common-index cohort\nn = {row.complete_indices_n:,} patients']
    reasons = [
        ('Subsequent ICU stays', row.excluded_nonfirst),
        ('Unavailable\nsource age', row.excluded_source_covariates),
        ('Age <18 years', row.excluded_age),
        ('ICU stay or hospital\nfollow-up <24 h', row.excluded_duration),
        ('Missing or invalid\nindex components\nin −6 to +24 h', row.excluded_laboratory)]
    counts = [int(row[c]) for c in ['all_icu_stays', 'first_patient_stays',
              'source_covariates_available', 'adult_first_stays',
              'eligible_pre_laboratory', 'complete_indices_n']]
    for i, (y, text) in enumerate(zip(ys, labels)):
        text_boxes.append(box(left, y, 27, 15, text,
                              fill if i == 5 else 'white', i == 5))
        if i == 5:
            continue
        reason, n = reasons[i]
        assert counts[i] - int(n) == counts[i+1]
        branch_y = (y+ys[i+1])/2
        arrow((left+13.5, y-7.5), (left+13.5, ys[i+1]+7.5))
        arrow((left+13.5, branch_y), (left+31, branch_y))
        text_boxes.append(box(left+31, branch_y, 26, 15,
                              f'{reason}\nn = {int(n):,}', fill, size=6.4))
        steps.append({'cohort': row.cohort, 'step': i+1,
                      'input_n': counts[i], 'excluded_n': int(n),
                      'retained_n': counts[i+1],
                      'exclusion_rule': reason.replace('\n', ' ')})

# Check actual text extents inside every box before exporting.
fig.canvas.draw()
renderer = fig.canvas.get_renderer()
for patch, label in text_boxes:
    pb = patch.get_window_extent(renderer)
    tb = label.get_window_extent(renderer)
    assert pb.contains(tb.x0, tb.y0) and pb.contains(tb.x1, tb.y1), label.get_text()
pd.DataFrame(steps).to_csv(S/'Figure1_patient_flow_steps.csv', index=False)
for ext in ['pdf', 'svg', 'png']:
    fig.savefig(F/('Figure1_Cohort_flowchart.'+ext), dpi=300)
plt.close(fig)
print('Figure 1: 3 cohorts, 15 sequential exclusions, all counts reconcile.')
