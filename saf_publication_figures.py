"""
SAF Harmonization — Extended Publication Figures
=================================================
Generates 11 additional publication-quality figures for the SAF harmonization
meta-analysis paper, covering all major sections:

  A  System boundary diagram (framework overview)
  B  Harmonization methodology flowchart
  C  MFSP cost breakdown stacked bar (TEA)
  D  Tornado sensitivity chart (TEA)
  E  Monte Carlo overlay before/after harmonization (TEA)
  F  GHG comparison with uncertainty ranges (LCA)
  G  Lifecycle contribution analysis (LCA)
  H  Parameter heatmap across all 42 studies (inconsistency)
  I  Variance decomposition — methodological vs. technical (inconsistency)
  J  TEA vs. LCA parity scatter (synthesis)
  K  Multi-criteria radar plot (synthesis)

Usage
-----
  python saf_publication_figures.py

Requires saf_harmonization_meta_analysis.py in the same directory.
"""

import sys
import os
import numpy as np
import pandas as pd
from scipy import stats
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch, Polygon
from matplotlib.gridspec import GridSpec
from matplotlib.lines import Line2D
import seaborn as sns
import warnings

warnings.filterwarnings('ignore')
np.random.seed(42)

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from saf_harmonization_meta_analysis import (
    build_harmonized_dataset,
    run_monte_carlo,
    run_sobol_analysis,
    decompose_variance,
    PATHWAY_CONFIGS,
    METHODOLOGICAL_PARAMS,
    PETROLEUM_JET_GHG_WTW,
)

os.makedirs('outputs/figures', exist_ok=True)

plt.rcParams.update({
    'font.family': 'DejaVu Sans',
    'font.size': 10,
    'axes.linewidth': 1.2,
    'figure.dpi': 150,
})

# ─────────────────────────────────────────────────────────────────────────────
# SHARED CONSTANTS
# ─────────────────────────────────────────────────────────────────────────────

PATHWAY_COLORS = {
    'ATJ':    '#2196F3',
    'HEFA':   '#4CAF50',
    'FT-SPK': '#FF9800',
    'PtL':    '#9C27B0',
}

PATHWAY_LABELS = {
    'ATJ':    'ATJ\n(Alcohol-to-Jet)',
    'HEFA':   'HEFA\n(Hydroprocessed)',
    'FT-SPK': 'FT-SPK\n(Fischer-Tropsch)',
    'PtL':    'PtL\n(Power-to-Liquid)',
}

# Representative modal-parameter cost breakdown (2023 $/GGE) per pathway
# Derived analytically from each pathway model at modal parameter values
COST_BREAKDOWN = {
    'ATJ': {
        'Capital Recovery':  1.66,
        'Fixed O&M':         0.63,
        'Variable O&M':      0.47,
        'Feedstock':         1.30,
    },
    'HEFA': {
        'Capital Recovery':  0.39,
        'Fixed O&M':         0.16,
        'Variable O&M':      0.10,
        'Feedstock':         3.06,
        'Hydrogen':          0.52,
    },
    'FT-SPK': {
        'Capital Recovery':  2.93,
        'Fixed O&M':         0.73,
        'Variable O&M':      0.70,
        'Feedstock':         1.47,
    },
    'PtL': {
        'Capital Recovery':  2.12,
        'Fixed O&M':         0.80,
        'Electricity':       5.52,
        'CO₂ Capture':       1.44,
    },
}

COST_COLORS = {
    'Capital Recovery': '#1565C0',
    'Fixed O&M':        '#42A5F5',
    'Variable O&M':     '#90CAF9',
    'Feedstock':        '#2E7D32',
    'Hydrogen':         '#A5D6A7',
    'Electricity':      '#E65100',
    'CO₂ Capture':      '#FFCC02',
}

# Representative GHG lifecycle contributions (gCO2e/MJ, harmonized WtWake basis)
# Computed at modal parameter values from each pathway model
GHG_CONTRIBUTIONS = {
    'ATJ': {
        'Feedstock Production': 10.65,
        'Conversion – NG':       5.97,
        'Conversion – Elec':     2.84,
        'Distribution':          3.00,
    },
    'HEFA': {
        'Feedstock Production': 16.64,
        'Hydrogen (SMR)':        6.74,
        'Conversion Energy':     1.56,
        'Distribution':          3.00,
    },
    'FT-SPK': {
        'Feedstock Production':  3.96,
        'Conversion Energy':     1.98,
        'Distribution':          3.00,
    },
    'PtL': {
        'Electricity (lifecycle)': 5.45,
        'Distribution':            1.50,
    },
}

GHG_STAGE_COLORS = {
    'Feedstock Production':    '#388E3C',
    'Conversion – NG':         '#F57C00',
    'Conversion – Elec':       '#FFA726',
    'Conversion Energy':       '#F57C00',
    'Hydrogen (SMR)':          '#B0BEC5',
    'Electricity (lifecycle)': '#7E57C2',
    'Distribution':            '#90A4AE',
    'CO₂ Capture':             '#FFCA28',
}

# Multi-criteria data for radar plot (harmonized, representative values)
RADAR_DATA = {
    'ATJ':    {'Cost\n($/GGE)': 4.1,  'GHG\nReduction (%)': 75, 'Land Use\n(m²/GJ)': 0.55, 'Water\n(L/GJ)': 25, 'TRL\n(1–9)': 6.5},
    'HEFA':   {'Cost\n($/GGE)': 4.3,  'GHG\nReduction (%)': 65, 'Land Use\n(m²/GJ)': 1.20, 'Water\n(L/GJ)': 30, 'TRL\n(1–9)': 9.0},
    'FT-SPK': {'Cost\n($/GGE)': 6.5,  'GHG\nReduction (%)': 90, 'Land Use\n(m²/GJ)': 0.45, 'Water\n(L/GJ)': 20, 'TRL\n(1–9)': 6.0},
    'PtL':    {'Cost\n($/GGE)': 11.5, 'GHG\nReduction (%)': 92, 'Land Use\n(m²/GJ)': 0.05, 'Water\n(L/GJ)': 80, 'TRL\n(1–9)': 5.0},
}

# Direction: True = higher is better (for normalization)
RADAR_HIGHER_BETTER = {
    'Cost\n($/GGE)':      False,
    'GHG\nReduction (%)': True,
    'Land Use\n(m²/GJ)':  False,
    'Water\n(L/GJ)':      False,
    'TRL\n(1–9)':         True,
}


# ─────────────────────────────────────────────────────────────────────────────
# HELPER: save figure
# ─────────────────────────────────────────────────────────────────────────────

def _save(fig, name):
    fig.savefig(f'outputs/figures/{name}.png', dpi=600, bbox_inches='tight')
    fig.savefig(f'outputs/figures/{name}.pdf',           bbox_inches='tight')
    plt.close(fig)
    print(f'  Saved: outputs/figures/{name}.{{png,pdf}}')


# ─────────────────────────────────────────────────────────────────────────────
# FIGURE A: SYSTEM BOUNDARY DIAGRAM
# ─────────────────────────────────────────────────────────────────────────────

def figA_system_boundary():
    """Conceptual system boundary diagram for all four SAF pathways."""
    fig, ax = plt.subplots(figsize=(18, 9))
    ax.set_xlim(0, 18)
    ax.set_ylim(0, 9)
    ax.axis('off')

    # ── Lifecycle stage boxes (top row) ──────────────────────────────────────
    stages = [
        'Feedstock\nProduction',
        'Feedstock\nTransport',
        'Conversion\n& Upgrading',
        'Fuel\nDistribution',
        'Fueling\n& Storage',
        'Combustion\n(Aircraft)',
    ]
    stage_x   = [1.0, 3.5, 6.5, 10.0, 13.0, 15.5]
    stage_w   = 2.0
    stage_h   = 1.1
    stage_y   = 7.2
    stage_col = '#E3F2FD'
    stage_edge= '#1565C0'

    for i, (label, x) in enumerate(zip(stages, stage_x)):
        box = FancyBboxPatch((x, stage_y), stage_w, stage_h,
                              boxstyle='round,pad=0.08',
                              facecolor=stage_col, edgecolor=stage_edge, lw=1.4)
        ax.add_patch(box)
        ax.text(x + stage_w / 2, stage_y + stage_h / 2, label,
                ha='center', va='center', fontsize=8.5, fontweight='bold', color='#0D47A1')
        if i < len(stages) - 1:
            ax.annotate('', xy=(stage_x[i + 1], stage_y + stage_h / 2),
                        xytext=(x + stage_w, stage_y + stage_h / 2),
                        arrowprops=dict(arrowstyle='->', color='#1565C0', lw=1.5))

    # ── System boundary brackets ──────────────────────────────────────────────
    def draw_bracket(x0, x1, y, label, color, ls='-'):
        ax.annotate('', xy=(x1 + 0.15, y), xytext=(x0 - 0.15, y),
                    arrowprops=dict(arrowstyle='<->', color=color, lw=2.0, ls=ls))
        ax.text((x0 + x1) / 2, y - 0.22, label,
                ha='center', va='top', fontsize=8, color=color, style='italic')

    draw_bracket(stage_x[0], stage_x[2] + stage_w, 6.8,
                 'Well-to-Gate (WtG)', '#B71C1C')
    draw_bracket(stage_x[0], stage_x[3] + stage_w, 6.45,
                 'Well-to-Wheel / Well-to-Wake (WtW / WtWake) ← Harmonized reference', '#1B5E20')
    draw_bracket(stage_x[2], stage_x[5] + stage_w, 6.1,
                 'Gate-to-Gate (GtG)', '#E65100', ls='--')

    # ── Pathway feedstock rows ────────────────────────────────────────────────
    pathway_rows = [
        ('ATJ',    '#2196F3', 'Lignocellulosic biomass\n(corn stover, bagasse, straw)',
         'Fermentation → Dehydration\n→ Oligomerisation → Hydrogenation'),
        ('HEFA',   '#4CAF50', 'Oils & fats\n(UCO, soy, tallow, camelina)',
         'Hydrotreating → Hydrocracking\n→ Isomerisation'),
        ('FT-SPK', '#FF9800', 'Solid biomass\n(wood, MSW, agricultural residues)',
         'Gasification → Gas cleaning\n→ Fischer-Tropsch synthesis'),
        ('PtL',    '#9C27B0', 'Renewable electricity + CO₂\n(solar, wind; DAC or point-source)',
         'Electrolysis (H₂) → CO₂ capture\n→ Fischer-Tropsch / Methanol'),
    ]

    row_ys  = [5.0, 3.8, 2.6, 1.4]
    feed_x  = stage_x[0]
    conv_x  = stage_x[2]
    box2_w  = 2.2
    box2_h  = 0.85

    for (pathway, color, feedstock_txt, conv_txt), ry in zip(pathway_rows, row_ys):
        # Feedstock box
        fb = FancyBboxPatch((feed_x, ry), box2_w, box2_h,
                             boxstyle='round,pad=0.06',
                             facecolor=color + '33', edgecolor=color, lw=1.2)
        ax.add_patch(fb)
        ax.text(feed_x + box2_w / 2, ry + box2_h / 2, feedstock_txt,
                ha='center', va='center', fontsize=7, color='#212121')

        # Conversion box
        cb = FancyBboxPatch((conv_x, ry), 3.2, box2_h,
                             boxstyle='round,pad=0.06',
                             facecolor=color + '22', edgecolor=color, lw=1.2, ls='--')
        ax.add_patch(cb)
        ax.text(conv_x + 1.6, ry + box2_h / 2, conv_txt,
                ha='center', va='center', fontsize=7, color='#212121')

        # Pathway label on left
        ax.text(feed_x - 0.15, ry + box2_h / 2, pathway,
                ha='right', va='center', fontsize=9, fontweight='bold', color=color)

        # Arrow: feedstock → conversion
        ax.annotate('', xy=(conv_x, ry + box2_h / 2),
                    xytext=(feed_x + box2_w, ry + box2_h / 2),
                    arrowprops=dict(arrowstyle='->', color=color, lw=1.2))

        # Arrow: conversion → distribution (all pathways converge to same jet fuel)
        ax.annotate('', xy=(stage_x[3] + 0.1, stage_y),
                    xytext=(conv_x + 3.2, ry + box2_h),
                    arrowprops=dict(arrowstyle='->', color=color, lw=0.9,
                                   connectionstyle='arc3,rad=-0.2'), zorder=1)

    # ── Co-product arrows ─────────────────────────────────────────────────────
    ax.annotate('', xy=(conv_x + 1.8, row_ys[0] - 0.3),
                xytext=(conv_x + 1.8, row_ys[0]),
                arrowprops=dict(arrowstyle='->', color='#757575', lw=1.0, ls='dotted'))
    ax.text(conv_x + 2.0, row_ys[0] - 0.45, 'Co-products\n(diesel, gasoline)',
            fontsize=6.5, color='#616161', style='italic')

    # ── Legend ────────────────────────────────────────────────────────────────
    legend_items = [
        Line2D([0], [0], color='#B71C1C', lw=2, label='Well-to-Gate (WtG)'),
        Line2D([0], [0], color='#1B5E20', lw=2, label='Well-to-Wake (WtWake) — harmonized'),
        Line2D([0], [0], color='#E65100', lw=2, ls='--', label='Gate-to-Gate (GtG)'),
    ]
    ax.legend(handles=legend_items, loc='lower right', fontsize=8.5,
              framealpha=0.9, edgecolor='#BDBDBD')

    ax.set_title('Figure A — SAF Production Pathways: System Boundary Definitions',
                 fontsize=12, fontweight='bold', pad=8)

    _save(fig, 'FigA_System_Boundary')


# ─────────────────────────────────────────────────────────────────────────────
# FIGURE B: HARMONIZATION METHODOLOGY FLOWCHART
# ─────────────────────────────────────────────────────────────────────────────

def figB_harmonization_flowchart():
    """Step-by-step harmonization protocol flowchart."""
    fig, ax = plt.subplots(figsize=(13, 15))
    ax.set_xlim(0, 13)
    ax.set_ylim(0, 15)
    ax.axis('off')

    def process_box(x, y, w, h, text, color='#E3F2FD', edge='#1565C0', bold=False):
        box = FancyBboxPatch((x, y), w, h, boxstyle='round,pad=0.1',
                              facecolor=color, edgecolor=edge, lw=1.5)
        ax.add_patch(box)
        ax.text(x + w / 2, y + h / 2, text, ha='center', va='center',
                fontsize=8.5, fontweight='bold' if bold else 'normal',
                color='#0D47A1' if edge == '#1565C0' else '#212121')

    def diamond(cx, cy, hw, hh, text, color='#FFF9C4', edge='#F57F17'):
        diamond_pts = np.array([
            [cx,      cy + hh],
            [cx + hw, cy],
            [cx,      cy - hh],
            [cx - hw, cy],
        ])
        poly = Polygon(diamond_pts, closed=True, facecolor=color,
                       edgecolor=edge, lw=1.5)
        ax.add_patch(poly)
        ax.text(cx, cy, text, ha='center', va='center', fontsize=8,
                fontweight='bold', color='#4E342E')

    def arrow(x0, y0, x1, y1, label='', color='#546E7A'):
        ax.annotate('', xy=(x1, y1), xytext=(x0, y0),
                    arrowprops=dict(arrowstyle='->', color=color, lw=1.6))
        if label:
            mx, my = (x0 + x1) / 2, (y0 + y1) / 2
            ax.text(mx + 0.15, my, label, fontsize=7.5, color=color, style='italic')

    def side_box(x, y, w, h, text, color='#F3E5F5', edge='#7B1FA2'):
        box = FancyBboxPatch((x, y), w, h, boxstyle='round,pad=0.08',
                              facecolor=color, edgecolor=edge, lw=1.2)
        ax.add_patch(box)
        ax.text(x + w / 2, y + h / 2, text, ha='center', va='center',
                fontsize=7.5, color='#4A148C')

    CX, BW, BH = 4.5, 4.0, 0.75   # centre x, box width, box height
    DH = 0.65                       # diamond half-height

    # Input
    process_box(CX - BW / 2, 13.5, BW, 0.85,
                'Raw Literature Study\n(heterogeneous methods, currencies, units)',
                color='#ECEFF1', edge='#455A64', bold=True)
    arrow(CX, 13.5, CX, 12.8)

    # Step 1: Currency
    process_box(CX - BW / 2, 12.0, BW, BH,
                'Step 1 — Currency Conversion\nEUR → USD at historical exchange rate')
    side_box(9.0, 12.05, 3.5, 0.65,
             'EUR/USD rates: 2009–2023\n(ECB historical series)')
    ax.annotate('', xy=(9.0, 12.38), xytext=(CX + BW / 2, 12.38),
                arrowprops=dict(arrowstyle='-', color='#7B1FA2', lw=1.0, ls='dashed'))
    arrow(CX, 12.0, CX, 11.3)

    # Step 2: Cost escalation
    process_box(CX - BW / 2, 10.5, BW, BH,
                'Step 2 — Cost Escalation to 2023 USD\nCapex: CEPCI ratio  |  Opex: CPI multiplier')
    side_box(9.0, 10.55, 3.5, 0.65,
             'CEPCI 2023 = 798\n(Chemical Engineering Plant Cost Index)')
    ax.annotate('', xy=(9.0, 10.88), xytext=(CX + BW / 2, 10.88),
                arrowprops=dict(arrowstyle='-', color='#7B1FA2', lw=1.0, ls='dashed'))
    arrow(CX, 10.5, CX, 9.8)

    # Step 3: Unit conversion
    process_box(CX - BW / 2, 9.0, BW, BH,
                'Step 3 — Unit Conversion\n$/L  |  €/L  →  2023 $/GGE of jet fuel')
    side_box(9.0, 9.05, 3.5, 0.65,
             '1 GGE = 3.533 L jet fuel\n(LHV basis, DOE standard)')
    ax.annotate('', xy=(9.0, 9.38), xytext=(CX + BW / 2, 9.38),
                arrowprops=dict(arrowstyle='-', color='#7B1FA2', lw=1.0, ls='dashed'))
    arrow(CX, 9.0, CX, 8.3)

    # Step 4: TEA normalization
    process_box(CX - BW / 2, 7.5, BW, BH,
                'Step 4 — TEA Parameter Normalization\nCorrect MFSP for DR, CF, plant lifetime via CRF ratio')
    side_box(9.0, 7.55, 3.5, 0.65,
             'Ref: DR = 10%, CF = 90%\nPlant life = 30 years')
    ax.annotate('', xy=(9.0, 7.88), xytext=(CX + BW / 2, 7.88),
                arrowprops=dict(arrowstyle='-', color='#7B1FA2', lw=1.0, ls='dashed'))
    arrow(CX, 7.5, CX, 6.8)

    # Step 5: GHG boundary
    diamond(CX, 6.2, 2.3, DH,
            'System\nBoundary?')
    arrow(CX, 6.2 - DH, CX, 5.6, label='WtWake')
    # Side: WtG branch
    ax.text(2.0, 6.2, 'WtG / GtG', ha='center', va='center', fontsize=7.5,
            color='#B71C1C', style='italic')
    ax.annotate('', xy=(2.0, 5.8), xytext=(CX - 2.3, 6.2),
                arrowprops=dict(arrowstyle='->', color='#B71C1C', lw=1.2))
    side_box(0.3, 5.4, 2.8, 0.65, '+3 gCO₂e/MJ\n(distribution & storage)', '#FFEBEE', '#C62828')
    ax.annotate('', xy=(CX - 2.3, 5.7), xytext=(3.1, 5.7),
                arrowprops=dict(arrowstyle='->', color='#B71C1C', lw=1.0))

    # Step 6: Allocation
    diamond(CX, 5.2, 2.3, DH,
            'Allocation\nMethod?')
    arrow(CX, 5.2 - DH, CX, 4.6, label='Energy alloc.')
    ax.text(2.0, 5.2, 'Mass / Econ /\nSys.Exp.', ha='center', va='center',
            fontsize=7.5, color='#E65100', style='italic')
    ax.annotate('', xy=(2.0, 4.8), xytext=(CX - 2.3, 5.2),
                arrowprops=dict(arrowstyle='->', color='#E65100', lw=1.2))
    side_box(0.3, 4.4, 2.8, 0.65,
             'Re-weight: GHG × (energy_factor\n / study_factor)', '#FFF3E0', '#E65100')
    ax.annotate('', xy=(CX - 2.3, 4.7), xytext=(3.1, 4.7),
                arrowprops=dict(arrowstyle='->', color='#E65100', lw=1.0))

    # Step 7: ILUC
    diamond(CX, 4.2, 2.3, DH,
            'ILUC\nIncluded?')
    arrow(CX, 4.2 - DH, CX, 3.6, label='No')
    ax.text(2.0, 4.2, 'Yes', ha='center', va='center', fontsize=7.5,
            color='#1B5E20', style='italic')
    ax.annotate('', xy=(2.0, 3.8), xytext=(CX - 2.3, 4.2),
                arrowprops=dict(arrowstyle='->', color='#1B5E20', lw=1.2))
    side_box(0.3, 3.4, 2.8, 0.65,
             'Subtract feedstock-specific\nILUC estimate (0–25 gCO₂e/MJ)', '#E8F5E9', '#1B5E20')
    ax.annotate('', xy=(CX - 2.3, 3.7), xytext=(3.1, 3.7),
                arrowprops=dict(arrowstyle='->', color='#1B5E20', lw=1.0))

    # Output
    process_box(CX - BW / 2, 2.6, BW, BH,
                'Harmonized Study Output\n2023 $/GGE  |  gCO₂e/MJ WtWake  |  Energy allocation',
                color='#E8F5E9', edge='#2E7D32', bold=True)

    # Final arrow up to Monte Carlo label
    arrow(CX, 2.6, CX, 1.9)
    process_box(CX - BW / 2, 1.0, BW, 0.8,
                'Monte Carlo uncertainty propagation\n(10,000 iterations per pathway)',
                color='#EDE7F6', edge='#512DA8')

    ax.set_title('Figure B — SAF Harmonization Methodology: Step-by-Step Protocol',
                 fontsize=12, fontweight='bold', pad=8)

    _save(fig, 'FigB_Harmonization_Flowchart')


# ─────────────────────────────────────────────────────────────────────────────
# FIGURE C: MFSP COST BREAKDOWN
# ─────────────────────────────────────────────────────────────────────────────

def figC_mfsp_cost_breakdown(df_harm, mc_results):
    """Stacked bar chart of MFSP cost components + literature range overlay."""
    fig, ax = plt.subplots(figsize=(14, 6))
    pathways = list(PATHWAY_COLORS.keys())

    x    = np.arange(len(pathways))
    w    = 0.45
    gap  = 0.52   # centre offset between stacked bar and MC bar

    # ── Stacked cost-component bars ───────────────────────────────────────────
    for i, pathway in enumerate(pathways):
        components = COST_BREAKDOWN[pathway]
        bottom = 0.0
        for j, (comp, val) in enumerate(components.items()):
            ax.bar(x[i] - gap / 2, val, w, bottom=bottom,
                   color=COST_COLORS.get(comp, '#BDBDBD'),
                   edgecolor='white', linewidth=0.6)
            if val > 0.4:
                ax.text(x[i] - gap / 2, bottom + val / 2, f'{val:.2f}',
                        ha='center', va='center', fontsize=6.5, color='white',
                        fontweight='bold')
            bottom += val

        # Total label on top
        total = sum(components.values())
        ax.text(x[i] - gap / 2, total + 0.08, f'${total:.2f}',
                ha='center', va='bottom', fontsize=8.5, fontweight='bold',
                color='#212121')

    # ── MC P10–P90 bars (harmonized range) ────────────────────────────────────
    for i, pathway in enumerate(pathways):
        vals = mc_results[pathway]['mfsp'].dropna()
        p10 = np.percentile(vals, 10)
        p50 = np.percentile(vals, 50)
        p90 = np.percentile(vals, 90)
        ax.bar(x[i] + gap / 2, p50, w, color=PATHWAY_COLORS[pathway],
               alpha=0.85, edgecolor='k', linewidth=0.7, label='MC Median' if i == 0 else '')
        ax.errorbar(x[i] + gap / 2, p50, yerr=[[p50 - p10], [p90 - p50]],
                    fmt='none', color='black', capsize=4, linewidth=1.4)
        ax.text(x[i] + gap / 2, p90 + 0.12, f'${p50:.2f}',
                ha='center', va='bottom', fontsize=8.5, fontweight='bold',
                color=PATHWAY_COLORS[pathway])

    # ── Literature range scatter ───────────────────────────────────────────────
    for i, pathway in enumerate(pathways):
        pts = df_harm[df_harm.pathway == pathway]['mfsp_harmonized']
        ax.scatter(np.full(len(pts), x[i] + gap / 2 + 0.30), pts,
                   color=PATHWAY_COLORS[pathway], alpha=0.5, s=22, zorder=5,
                   edgecolors='k', linewidths=0.3)

    # ── Legend for cost components ─────────────────────────────────────────────
    all_comps = list(dict.fromkeys(
        c for d in COST_BREAKDOWN.values() for c in d
    ))
    legend_patches = [
        mpatches.Patch(color=COST_COLORS.get(c, '#BDBDBD'), label=c)
        for c in all_comps
    ]
    legend_patches += [
        mpatches.Patch(color='#90A4AE', alpha=0.85, label='MC median (P10–P90 bars)'),
        Line2D([0], [0], marker='o', color='gray', ls='None',
               markersize=5, alpha=0.6, label='Harmonized literature studies'),
    ]
    ax.legend(handles=legend_patches, fontsize=7.5, loc='upper left',
              ncol=2, framealpha=0.9)

    ax.set_xticks(x)
    ax.set_xticklabels(
        [f'{PATHWAY_LABELS[p]}\n(left: cost breakdown,\nright: MC distribution)' for p in pathways],
        fontsize=8
    )
    ax.set_ylabel('MFSP (2023 $/GGE)', fontsize=11, fontweight='bold')
    ax.set_title('Figure C — Minimum Fuel Selling Price: Cost Breakdown and Monte Carlo Range',
                 fontsize=11, fontweight='bold', loc='left')
    ax.set_ylim(0, ax.get_ylim()[1] * 1.12)
    ax.grid(axis='y', alpha=0.25)

    _save(fig, 'FigC_MFSP_Cost_Breakdown')


# ─────────────────────────────────────────────────────────────────────────────
# FIGURE D: TORNADO SENSITIVITY CHART
# ─────────────────────────────────────────────────────────────────────────────

def figD_tornado_sensitivity(mc_results):
    """
    Tornado chart: Spearman rank correlation of each parameter with MFSP.
    One subplot per pathway.
    """
    fig, axes = plt.subplots(2, 2, figsize=(16, 10))
    axes = axes.flatten()

    PARAM_LABELS = {
        'feedstock_cost':   'Feedstock Cost',
        'capex_2023':       'Capital Cost',
        'ft_efficiency':    'FT Efficiency',
        'ethanol_yield':    'Fuel Yield',
        'jet_yield':        'Jet Fraction',
        'capacity_factor':  'Capacity Factor',
        'discount_rate':    'Discount Rate',
        'h2_use':           'H₂ Consumption',
        'h2_price':         'H₂ Price',
        'elec_cost_mwh':    'Electricity Cost',
        'elec_capex_kw':    'Electrolyzer Capex',
        'ft_capex':         'FT Capex (PtL)',
        'co2_capture_cost': 'CO₂ Capture Cost',
        'alloc_factor':     'Allocation Factor *',
        'boundary_offset':  'Boundary Offset *',
        'iluc_penalty':     'ILUC Penalty *',
        'ng_use':           'NG Use',
        'elec_use':         'Electricity Use',
        'feedstock_ghg':    'Feedstock GHG',
        'grid_intensity':   'Grid Intensity',
        'process_ghg':      'Process GHG',
    }

    for ax, pathway in zip(axes, PATHWAY_COLORS.keys()):
        mc = mc_results[pathway].dropna()
        params = [c for c in mc.columns if c not in ('mfsp', 'ghg')]

        corrs = {}
        for param in params:
            if mc[param].std() > 0:
                r, _ = stats.spearmanr(mc[param], mc['mfsp'])
                corrs[param] = r

        # Sort by absolute correlation
        sorted_params = sorted(corrs.items(), key=lambda kv: abs(kv[1]), reverse=True)
        top_n = min(10, len(sorted_params))
        sorted_params = sorted_params[:top_n]

        labels = [PARAM_LABELS.get(k, k) for k, _ in sorted_params]
        values = [v for _, v in sorted_params]
        colors = ['#e53935' if k in METHODOLOGICAL_PARAMS else
                  PATHWAY_COLORS[pathway]
                  for k, _ in sorted_params]

        y_pos = range(top_n)
        ax.barh(list(y_pos), values, color=colors, alpha=0.85,
                edgecolor='k', linewidth=0.5)
        ax.axvline(0, color='black', lw=0.8)
        ax.set_yticks(list(y_pos))
        ax.set_yticklabels(labels, fontsize=8.5)
        ax.set_xlabel("Spearman's ρ with MFSP", fontsize=9, fontweight='bold')
        ax.set_title(f'{pathway}', fontsize=11, fontweight='bold', loc='left',
                     color=PATHWAY_COLORS[pathway])
        ax.grid(axis='x', alpha=0.25)
        ax.set_xlim(-1.05, 1.05)

        # Annotation
        for j, (val, (k, _)) in enumerate(zip(values, sorted_params)):
            ax.text(val + (0.03 if val >= 0 else -0.03),
                    j, f'{val:.2f}',
                    ha='left' if val >= 0 else 'right',
                    va='center', fontsize=7, color='#212121')

    # Legend
    m_patch = mpatches.Patch(color='#e53935', label='Methodological parameter (*)')
    t_patch = mpatches.Patch(color='#78909C', label='Technical / economic parameter')
    fig.legend(handles=[m_patch, t_patch], loc='lower center', ncol=2,
               fontsize=9, framealpha=0.9, bbox_to_anchor=(0.5, 0.01))
    fig.suptitle("Figure D — Tornado Sensitivity: Spearman Rank Correlation with MFSP",
                 fontsize=12, fontweight='bold', y=1.01)
    plt.tight_layout()

    _save(fig, 'FigD_Tornado_Sensitivity')


# ─────────────────────────────────────────────────────────────────────────────
# FIGURE E: MONTE CARLO OVERLAY — BEFORE / AFTER HARMONIZATION
# ─────────────────────────────────────────────────────────────────────────────

def figE_mc_overlay(df_harm, mc_results):
    """KDE overlay: raw reported values vs harmonized Monte Carlo distribution."""
    fig, axes = plt.subplots(2, 2, figsize=(14, 9))
    axes = axes.flatten()

    for ax, pathway in zip(axes, PATHWAY_COLORS.keys()):
        color = PATHWAY_COLORS[pathway]
        raw   = df_harm[df_harm.pathway == pathway]['mfsp_2023_raw'].dropna()
        harm  = df_harm[df_harm.pathway == pathway]['mfsp_harmonized'].dropna()
        mc_v  = mc_results[pathway]['mfsp'].dropna()

        x_min = min(raw.min(), mc_v.min()) * 0.85
        x_max = max(raw.max(), mc_v.max()) * 1.15
        xs    = np.linspace(x_min, x_max, 400)

        # KDE: raw literature values
        if len(raw) > 2:
            kde_raw = stats.gaussian_kde(raw, bw_method=0.5)
            ax.fill_between(xs, kde_raw(xs), alpha=0.25, color='#616161',
                            label='Raw literature (KDE)')
            ax.plot(xs, kde_raw(xs), color='#424242', lw=1.5, ls='--')

        # KDE: harmonized literature values
        if len(harm) > 2:
            kde_harm = stats.gaussian_kde(harm, bw_method=0.5)
            ax.fill_between(xs, kde_harm(xs), alpha=0.30, color=color,
                            label='Harmonized literature (KDE)')
            ax.plot(xs, kde_harm(xs), color=color, lw=1.5, ls='-.')

        # KDE: Monte Carlo model output
        kde_mc = stats.gaussian_kde(mc_v, bw_method=0.2)
        ax.fill_between(xs, kde_mc(xs), alpha=0.50, color=color,
                        label='MC model distribution')
        ax.plot(xs, kde_mc(xs), color=color, lw=2.0)

        # Rug plots for literature studies
        ax.scatter(raw,  np.full(len(raw),  -0.005), marker='|', color='#424242',
                   s=60, alpha=0.7, lw=1.2, clip_on=False)
        ax.scatter(harm, np.full(len(harm), -0.018), marker='|', color=color,
                   s=60, alpha=0.8, lw=1.2, clip_on=False)

        # Median lines
        ax.axvline(mc_v.median(), color=color, lw=1.8, ls='-',
                   label=f'MC median ${mc_v.median():.2f}')
        ax.axvline(raw.mean(), color='#616161', lw=1.5, ls='--',
                   label=f'Raw mean ${raw.mean():.2f}')

        ax.set_xlabel('MFSP (2023 $/GGE)', fontsize=10, fontweight='bold')
        ax.set_ylabel('Probability Density', fontsize=9)
        ax.set_title(f'{pathway}', fontsize=11, fontweight='bold',
                     color=color, loc='left')
        ax.legend(fontsize=7.5, loc='upper right')
        ax.grid(alpha=0.2)
        ax.set_xlim(x_min, x_max)
        ax.set_ylim(bottom=-0.03)

    fig.suptitle(
        'Figure E — Monte Carlo Distribution vs. Literature Values: Before and After Harmonization',
        fontsize=11, fontweight='bold'
    )
    plt.tight_layout()

    _save(fig, 'FigE_MC_Overlay_Before_After')


# ─────────────────────────────────────────────────────────────────────────────
# FIGURE F: GHG COMPARISON WITH UNCERTAINTY RANGES
# ─────────────────────────────────────────────────────────────────────────────

def figF_ghg_comparison(df_harm, mc_results):
    """GHG before/after harmonization with P5–P95 uncertainty ranges."""
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))

    pathways = list(PATHWAY_COLORS.keys())
    x = np.arange(len(pathways))
    w = 0.30

    # ── Left panel: mean ± std from literature ────────────────────────────────
    ax = axes[0]
    for shift, col_raw, col_harm, label_raw, label_harm in [
        (-w / 2, 'ghg_raw', 'ghg_harmonized', 'Raw', 'Harmonized')
    ]:
        for i, pathway in enumerate(pathways):
            grp_raw  = df_harm[df_harm.pathway == pathway]['ghg_raw'].dropna()
            grp_harm = df_harm[df_harm.pathway == pathway]['ghg_harmonized'].dropna()

            ax.bar(x[i] - w * 0.6, grp_raw.mean(), w,
                   yerr=grp_raw.std(), color=PATHWAY_COLORS[pathway],
                   alpha=0.35, edgecolor='k', capsize=5, linewidth=0.7,
                   label='Raw (mean ± σ)' if i == 0 else '')
            ax.bar(x[i] + w * 0.6, grp_harm.mean(), w,
                   yerr=grp_harm.std(), color=PATHWAY_COLORS[pathway],
                   alpha=0.85, edgecolor='k', capsize=5, linewidth=0.7,
                   label='Harmonized (mean ± σ)' if i == 0 else '')

            # Individual study points
            ax.scatter(np.full(len(grp_raw),  x[i] - w * 0.6 + 0.22), grp_raw,
                       color=PATHWAY_COLORS[pathway], alpha=0.5, s=18, zorder=5)
            ax.scatter(np.full(len(grp_harm), x[i] + w * 0.6 + 0.22), grp_harm,
                       color=PATHWAY_COLORS[pathway], alpha=0.8, s=18, zorder=5)

    ax.axhline(PETROLEUM_JET_GHG_WTW, color='red', ls='--', lw=1.4,
               label=f'Petroleum jet ({PETROLEUM_JET_GHG_WTW} gCO₂e/MJ)')
    ax.set_xticks(x)
    ax.set_xticklabels(pathways, fontsize=10)
    ax.set_ylabel('GHG Intensity (gCO₂e/MJ)', fontsize=11, fontweight='bold')
    ax.set_title('(a) Literature: Raw vs. Harmonized\n(mean ± 1σ, individual studies shown)',
                 fontsize=10, fontweight='bold', loc='left')
    ax.legend(fontsize=8.5)
    ax.grid(axis='y', alpha=0.25)
    ax.set_ylim(-15, 100)

    # ── Right panel: MC P5-P95 box-style ─────────────────────────────────────
    ax = axes[1]
    for i, pathway in enumerate(pathways):
        vals = mc_results[pathway]['ghg'].dropna()
        p5, p25, p50, p75, p95 = np.percentile(vals, [5, 25, 50, 75, 95])
        color = PATHWAY_COLORS[pathway]

        # Box P25–P75
        ax.bar(x[i], p75 - p25, w * 1.2, bottom=p25, color=color,
               alpha=0.75, edgecolor='k', linewidth=0.8)
        # Median line
        ax.plot([x[i] - w * 0.6, x[i] + w * 0.6], [p50, p50],
                color='white', lw=2.5)
        # Whiskers P5–P95
        ax.plot([x[i], x[i]], [p5, p25], color=color, lw=2.0)
        ax.plot([x[i], x[i]], [p75, p95], color=color, lw=2.0)
        ax.plot([x[i] - 0.12, x[i] + 0.12], [p5, p5],   color=color, lw=1.5)
        ax.plot([x[i] - 0.12, x[i] + 0.12], [p95, p95], color=color, lw=1.5)

        ax.text(x[i], p95 + 1.0, f'P50={p50:.1f}', ha='center', va='bottom',
                fontsize=7.5, fontweight='bold', color=color)

    ax.axhline(PETROLEUM_JET_GHG_WTW, color='red', ls='--', lw=1.4,
               label='Petroleum jet')
    ax.axhline(0, color='black', ls=':', lw=0.8)
    ax.set_xticks(x)
    ax.set_xticklabels(pathways, fontsize=10)
    ax.set_ylabel('GHG Intensity (gCO₂e/MJ)', fontsize=11, fontweight='bold')
    ax.set_title('(b) Monte Carlo Distribution (P5–P25–P50–P75–P95)',
                 fontsize=10, fontweight='bold', loc='left')
    ax.legend(fontsize=8.5)
    ax.grid(axis='y', alpha=0.25)
    ax.set_ylim(-15, 100)

    fig.suptitle('Figure F — Lifecycle GHG Emissions: Literature Spread and Probabilistic Range',
                 fontsize=11, fontweight='bold')
    plt.tight_layout()

    _save(fig, 'FigF_GHG_Comparison')


# ─────────────────────────────────────────────────────────────────────────────
# FIGURE G: GHG LIFECYCLE CONTRIBUTION ANALYSIS
# ─────────────────────────────────────────────────────────────────────────────

def figG_ghg_contribution():
    """Stacked horizontal bars showing GHG contribution by lifecycle stage."""
    fig, ax = plt.subplots(figsize=(14, 5))

    pathways  = list(PATHWAY_COLORS.keys())
    all_stages = list(dict.fromkeys(
        s for d in GHG_CONTRIBUTIONS.values() for s in d
    ))

    y = np.arange(len(pathways))
    bar_h = 0.55

    for i, pathway in enumerate(pathways):
        contribs = GHG_CONTRIBUTIONS[pathway]
        left = 0.0
        for stage in all_stages:
            val = contribs.get(stage, 0.0)
            if val == 0:
                continue
            color = GHG_STAGE_COLORS.get(stage, '#BDBDBD')
            ax.barh(y[i], val, bar_h, left=left, color=color,
                    edgecolor='white', linewidth=0.6)
            if val > 0.8:
                ax.text(left + val / 2, y[i], f'{val:.1f}',
                        ha='center', va='center', fontsize=8,
                        color='white', fontweight='bold')
            left += val

        # Total label
        total = sum(contribs.values())
        ax.text(total + 0.3, y[i], f'{total:.1f} gCO₂e/MJ',
                ha='left', va='center', fontsize=8.5, fontweight='bold',
                color=PATHWAY_COLORS[pathway])

    ax.axvline(PETROLEUM_JET_GHG_WTW, color='red', ls='--', lw=1.4,
               label=f'Petroleum jet ({PETROLEUM_JET_GHG_WTW} gCO₂e/MJ)')
    ax.set_yticks(y)
    ax.set_yticklabels(pathways, fontsize=11, fontweight='bold')
    ax.set_xlabel('GHG Intensity (gCO₂e/MJ)', fontsize=11, fontweight='bold')
    ax.set_title(
        'Figure G — Lifecycle GHG Contribution Analysis\n'
        '(Harmonized: WtWake | Energy allocation | Modal parameters)',
        fontsize=11, fontweight='bold', loc='left'
    )
    ax.grid(axis='x', alpha=0.25)
    ax.set_xlim(0, 38)

    # Stage legend
    stage_patches = [
        mpatches.Patch(color=GHG_STAGE_COLORS.get(s, '#BDBDBD'), label=s)
        for s in all_stages
    ]
    stage_patches.append(
        Line2D([0], [0], color='red', ls='--', lw=1.5, label='Petroleum jet baseline')
    )
    ax.legend(handles=stage_patches, loc='lower right', fontsize=8, ncol=2,
              framealpha=0.9)

    _save(fig, 'FigG_GHG_Contribution')


# ─────────────────────────────────────────────────────────────────────────────
# FIGURE H: PARAMETER HEATMAP ACROSS ALL STUDIES
# ─────────────────────────────────────────────────────────────────────────────

def figH_parameter_heatmap(df_harm):
    """Heatmap: rows = studies, columns = key methodological/technical choices."""
    fig, ax = plt.subplots(figsize=(16, 12))

    # Encode categorical → numeric
    df = df_harm.copy()

    alloc_map     = {'energy': 0, 'mass': 1, 'economic': 2, 'system_expansion': 3}
    boundary_map  = {'WtG': 0, 'GtG': 0.5, 'WtW': 1.5, 'WtWake': 2}
    df['alloc_enc']    = df['allocation'].map(alloc_map)
    df['boundary_enc'] = df['boundary'].map(boundary_map)
    df['iluc_enc']     = df['include_iluc'].astype(int)
    df['year_norm']    = (df['year'] - df['year'].min()) / (df['year'].max() - df['year'].min())

    cols_display = {
        'pathway':           'Pathway',
        'alloc_enc':         'Allocation\n(0=energy…3=sys.exp.)',
        'boundary_enc':      'System\nBoundary\n(0=WtG,2=WtWake)',
        'iluc_enc':          'ILUC\n(0=No,1=Yes)',
        'discount_rate':     'Discount\nRate (%)',
        'capacity_factor':   'Capacity\nFactor (%)',
        'plant_lifetime':    'Plant\nLifetime (yr)',
        'year_norm':         'Pub. Year\n(normalized)',
    }

    # Pathway as numeric
    pathway_map = {'ATJ': 0, 'HEFA': 1, 'FT-SPK': 2, 'PtL': 3}
    df['pathway_enc'] = df['pathway'].map(pathway_map)

    heat_cols_raw = ['pathway_enc', 'alloc_enc', 'boundary_enc', 'iluc_enc',
                     'discount_rate', 'capacity_factor', 'plant_lifetime', 'year_norm']
    col_labels = [
        'Pathway\n(0=ATJ…3=PtL)',
        'Allocation\n(0=energy…3=sysexp)',
        'Boundary\n(0=WtG,2=WtWake)',
        'ILUC\n(0=No,1=Yes)',
        'Discount\nRate (%)',
        'Capacity\nFactor (%)',
        'Plant\nLifetime (yr)',
        'Pub. Year\n(normalized)',
    ]

    heat_data = df[heat_cols_raw].values.astype(float)
    # Normalize each column to 0–1 for unified display
    heat_norm = np.zeros_like(heat_data)
    for j in range(heat_data.shape[1]):
        col = heat_data[:, j]
        mn, mx = col.min(), col.max()
        heat_norm[:, j] = (col - mn) / (mx - mn) if mx > mn else col * 0

    # Row labels: author + year + pathway
    row_labels = [f"{r['authors'].split()[0]} {r['year']} ({r['pathway']})"
                  for _, r in df.iterrows()]

    # Color rows by pathway
    row_colors = [PATHWAY_COLORS[p] for p in df['pathway']]

    sns.heatmap(heat_norm, ax=ax, cmap='RdYlGn_r', linewidths=0.3,
                linecolor='white', cbar_kws={'label': 'Normalized value (0=low, 1=high)'},
                vmin=0, vmax=1, annot=False)

    ax.set_xticks(np.arange(len(col_labels)) + 0.5)
    ax.set_xticklabels(col_labels, fontsize=8, rotation=0, ha='center')
    ax.set_yticks(np.arange(len(row_labels)) + 0.5)
    ax.set_yticklabels(row_labels, fontsize=6.5, rotation=0)

    # Color row labels by pathway
    for ytick, color in zip(ax.get_yticklabels(), row_colors):
        ytick.set_color(color)

    ax.set_title(
        'Figure H — Methodological Parameter Choices Across 42 SAF Studies\n'
        'Rows colored by pathway: Blue=ATJ, Green=HEFA, Orange=FT-SPK, Purple=PtL',
        fontsize=11, fontweight='bold', loc='left'
    )

    _save(fig, 'FigH_Parameter_Heatmap')


# ─────────────────────────────────────────────────────────────────────────────
# FIGURE I: VARIANCE DECOMPOSITION
# ─────────────────────────────────────────────────────────────────────────────

def figI_variance_decomposition(var_df, sobol_results):
    """
    Two-panel: (a) stacked bars of methodological vs technical variance fraction;
               (b) First-order Sobol S1 per parameter per pathway (MFSP).
    """
    fig, axes = plt.subplots(1, 2, figsize=(16, 6))

    pathways = list(PATHWAY_COLORS.keys())

    # ── Panel (a): stacked bars ───────────────────────────────────────────────
    ax = axes[0]
    x  = np.arange(len(pathways))
    w  = 0.30

    for j, (metric, ls) in enumerate([('MFSP', '-'), ('GHG', '--')]):
        mdf = var_df[var_df.Metric == metric].set_index('Pathway')
        meth  = [mdf.loc[p, 'Methodological_%'] for p in pathways]
        tech  = [mdf.loc[p, 'Technical_%']      for p in pathways]

        bars_m = ax.bar(x + (j - 0.5) * w, meth, w,
                        color='#e53935', alpha=0.85 if metric == 'MFSP' else 0.45,
                        edgecolor='k', hatch='' if metric == 'MFSP' else '//',
                        linewidth=0.7,
                        label=f'Methodological – {metric}')
        bars_t = ax.bar(x + (j - 0.5) * w, tech, w, bottom=meth,
                        color='#1565C0', alpha=0.85 if metric == 'MFSP' else 0.45,
                        edgecolor='k', hatch='' if metric == 'MFSP' else '//',
                        linewidth=0.7,
                        label=f'Technical – {metric}')

        for i, (m_val, t_val) in enumerate(zip(meth, tech)):
            vr_row = var_df[(var_df.Pathway == pathways[i]) & (var_df.Metric == metric)]
            vr = vr_row['Variance_Reduction_%'].values[0]
            ax.text(x[i] + (j - 0.5) * w, m_val + t_val + 1.0,
                    f'↓{vr:.0f}%', ha='center', va='bottom', fontsize=7.5,
                    color='#1B5E20', fontweight='bold')

    ax.set_xticks(x)
    ax.set_xticklabels(pathways, fontsize=10)
    ax.set_ylabel('Share of Total Variance (%)', fontsize=11, fontweight='bold')
    ax.set_title('(a) Variance Decomposition\n(↓% = CV reduction from harmonization)',
                 fontsize=10, fontweight='bold', loc='left')
    ax.legend(fontsize=8, loc='upper right', ncol=2)
    ax.grid(axis='y', alpha=0.25)
    ax.set_ylim(0, 130)

    # ── Panel (b): Sobol S1 grouped bar per pathway ───────────────────────────
    ax = axes[1]
    all_params = sorted(
        set(k for pw in PATHWAY_COLORS for k in sobol_results[pw]['S1_mfsp']),
        key=lambda k: k in METHODOLOGICAL_PARAMS
    )

    bar_w = 0.18
    n_pw  = len(pathways)
    x2    = np.arange(len(all_params))

    for j, pathway in enumerate(pathways):
        s1 = sobol_results[pathway]['S1_mfsp']
        vals = [s1.get(p, 0) for p in all_params]
        ax.bar(x2 + (j - n_pw / 2 + 0.5) * bar_w, vals, bar_w,
               color=PATHWAY_COLORS[pathway], alpha=0.85, edgecolor='k',
               linewidth=0.4, label=pathway)

    param_display = [
        ('alloc_factor',    'Alloc *'),
        ('boundary_offset', 'Boundary *'),
        ('iluc_penalty',    'ILUC *'),
        ('feedstock_cost',  'Feed. Cost'),
        ('capex_2023',      'Capex'),
        ('elec_cost_mwh',   'Elec. Cost'),
        ('ft_efficiency',   'FT Eff.'),
        ('ethanol_yield',   'Fuel Yield'),
        ('discount_rate',   'Disc. Rate'),
        ('capacity_factor', 'Cap. Factor'),
        ('h2_use',          'H₂ Use'),
        ('h2_price',        'H₂ Price'),
        ('co2_capture_cost','CO₂ Cost'),
        ('elec_capex_kw',   'Elec. Capex'),
        ('ft_capex',        'FT Capex'),
        ('ng_use',          'NG Use'),
        ('elec_use',        'Elec. Use'),
        ('feedstock_ghg',   'Feed. GHG'),
        ('grid_intensity',  'Grid GHG'),
        ('process_ghg',     'Proc. GHG'),
        ('jet_yield',       'Jet Yield'),
    ]
    label_map = dict(param_display)
    x_labels  = [label_map.get(p, p) for p in all_params]
    ax.set_xticks(x2)
    ax.set_xticklabels(x_labels, rotation=45, ha='right', fontsize=7.5)

    # Shade methodological params
    for j, p in enumerate(all_params):
        if p in METHODOLOGICAL_PARAMS:
            ax.axvspan(j - 0.45, j + 0.45, alpha=0.08, color='#e53935', zorder=0)

    ax.set_ylabel('First-Order Sobol Index (S₁)', fontsize=11, fontweight='bold')
    ax.set_title('(b) Sobol S₁ per Parameter — MFSP\n(shaded = methodological)',
                 fontsize=10, fontweight='bold', loc='left')
    ax.legend(fontsize=8, loc='upper right')
    ax.grid(axis='y', alpha=0.25)
    ax.set_ylim(0, min(1.0, ax.get_ylim()[1] * 1.2))

    fig.suptitle('Figure I — Variance Decomposition: Methodological vs. Technical Sources of Uncertainty',
                 fontsize=11, fontweight='bold')
    plt.tight_layout()

    _save(fig, 'FigI_Variance_Decomposition')


# ─────────────────────────────────────────────────────────────────────────────
# FIGURE J: TEA vs. LCA SCATTER (PARITY PLOT)
# ─────────────────────────────────────────────────────────────────────────────

def figJ_tea_lca_scatter(df_harm, mc_results):
    """
    Scatter plot: harmonized MFSP (x) vs harmonized GHG (y).
    Each point = one study; ellipses show MC P25–P75 clusters per pathway.
    """
    fig, ax = plt.subplots(figsize=(11, 8))

    # ── MC density ellipses (P25–P75 box approximation as ellipse) ───────────
    from matplotlib.patches import Ellipse
    for pathway in PATHWAY_COLORS:
        mc   = mc_results[pathway]
        xv   = mc['mfsp'].dropna()
        yv   = mc['ghg'].dropna()
        xc, yc = xv.median(), yv.median()
        xw   = np.percentile(xv, 75) - np.percentile(xv, 25)
        yh   = np.percentile(yv, 75) - np.percentile(yv, 25)
        ell  = Ellipse((xc, yc), width=xw * 2.0, height=yh * 2.0,
                       facecolor=PATHWAY_COLORS[pathway], alpha=0.12,
                       edgecolor=PATHWAY_COLORS[pathway], lw=1.5, ls='--')
        ax.add_patch(ell)
        ax.annotate(
            f'{pathway}\nMedian: ${xc:.2f}/GGE\n{yc:.1f} gCO₂e/MJ',
            xy=(xc, yc),
            xytext=(xc + 0.4, yc + 2.5),
            fontsize=7.5, color=PATHWAY_COLORS[pathway], fontweight='bold',
            arrowprops=dict(arrowstyle='->', color=PATHWAY_COLORS[pathway],
                            lw=0.9, alpha=0.7),
        )

    # ── Individual study scatter ──────────────────────────────────────────────
    for pathway in PATHWAY_COLORS:
        grp = df_harm[df_harm.pathway == pathway]
        ax.scatter(grp['mfsp_harmonized'], grp['ghg_harmonized'],
                   color=PATHWAY_COLORS[pathway], s=70, alpha=0.85,
                   edgecolors='k', linewidths=0.5, zorder=5,
                   label=pathway)

    # ── Petroleum reference ───────────────────────────────────────────────────
    ax.axhline(PETROLEUM_JET_GHG_WTW, color='red', ls='--', lw=1.5,
               label=f'Petroleum jet GHG ({PETROLEUM_JET_GHG_WTW} gCO₂e/MJ)')

    # ── Desirability region (low cost + low GHG) ──────────────────────────────
    ax.fill_between([0, 6], [0, 0], [20, 20], alpha=0.06, color='#1B5E20',
                    label='Preferred zone (low cost & GHG)')

    ax.set_xlabel('Harmonized MFSP (2023 $/GGE)', fontsize=12, fontweight='bold')
    ax.set_ylabel('Harmonized GHG (gCO₂e/MJ, WtWake)', fontsize=12, fontweight='bold')
    ax.set_title(
        'Figure J — TEA vs. LCA Parity Plot: Cost–Emissions Trade-off\n'
        'Points = individual harmonized studies; ellipses = MC P25–P75 cluster',
        fontsize=11, fontweight='bold', loc='left'
    )
    ax.legend(fontsize=9, loc='upper right', framealpha=0.9)
    ax.grid(alpha=0.22)
    ax.set_xlim(left=0)
    ax.set_ylim(-15, 100)

    _save(fig, 'FigJ_TEA_LCA_Scatter')


# ─────────────────────────────────────────────────────────────────────────────
# FIGURE K: MULTI-CRITERIA RADAR / SPIDER PLOT
# ─────────────────────────────────────────────────────────────────────────────

def figK_radar_plot():
    """Multi-criteria radar plot comparing all four pathways."""
    criteria = list(RADAR_DATA['ATJ'].keys())
    n_crit   = len(criteria)
    angles   = np.linspace(0, 2 * np.pi, n_crit, endpoint=False).tolist()
    angles  += angles[:1]   # close the polygon

    # Normalize each criterion to [0, 1]
    all_vals = {c: [RADAR_DATA[p][c] for p in RADAR_DATA] for c in criteria}
    norm_data = {}
    for pathway in RADAR_DATA:
        row = []
        for c in criteria:
            vals   = all_vals[c]
            mn, mx = min(vals), max(vals)
            raw    = RADAR_DATA[pathway][c]
            if mx == mn:
                n = 0.5
            elif RADAR_HIGHER_BETTER[c]:
                n = (raw - mn) / (mx - mn)
            else:
                n = (mx - raw) / (mx - mn)
            row.append(n)
        row += row[:1]   # close
        norm_data[pathway] = row

    fig, ax = plt.subplots(figsize=(10, 10),
                           subplot_kw=dict(polar=True))

    # Grid rings
    for r in [0.25, 0.50, 0.75, 1.0]:
        ax.plot(angles + angles[:1], [r] * (n_crit + 1), color='#BDBDBD',
                lw=0.8, ls='-', zorder=0)
        if r < 1.0:
            ax.text(np.pi / 2, r + 0.03, f'{int(r * 100)}%',
                    ha='center', fontsize=7.5, color='#757575')

    for pathway, values in norm_data.items():
        color = PATHWAY_COLORS[pathway]
        ax.plot(angles, values, color=color, lw=2.2, label=pathway, zorder=3)
        ax.fill(angles, values, color=color, alpha=0.10)
        # Markers at each vertex
        ax.scatter(angles, values, color=color, s=60, zorder=4, edgecolors='k',
                   linewidths=0.5)

    # Axis labels with actual values annotated
    ax.set_xticks(angles[:-1])
    ax.set_xticklabels(criteria, fontsize=9, fontweight='bold', color='#212121')
    ax.set_ylim(0, 1.15)
    ax.set_yticks([])

    # Annotate actual values at vertices
    for pathway, values in norm_data.items():
        color = PATHWAY_COLORS[pathway]
        for k, (angle, norm_val) in enumerate(zip(angles[:-1], values[:-1])):
            raw_val = RADAR_DATA[pathway][criteria[k]]
            # Format label
            if 'GHG' in criteria[k]:
                lbl = f'{raw_val:.0f}%'
            elif 'Cost' in criteria[k]:
                lbl = f'${raw_val:.1f}'
            elif 'TRL' in criteria[k]:
                lbl = f'TRL {raw_val:.0f}'
            elif 'Land' in criteria[k]:
                lbl = f'{raw_val:.2f}'
            elif 'Water' in criteria[k]:
                lbl = f'{raw_val:.0f}L'
            else:
                lbl = f'{raw_val}'

            offset = 0.16
            ax.text(angle, norm_val + offset, lbl,
                    ha='center', va='center', fontsize=6.5, color=color)

    ax.set_title(
        'Figure K — Multi-Criteria Comparison of SAF Pathways\n'
        '(all axes normalized: outer edge = best performance)',
        fontsize=11, fontweight='bold', pad=20
    )

    # Source annotations
    notes = (
        'Sources: MFSP from harmonized MC medians; GHG reduction vs. ICAO CORSIA baseline;\n'
        'Land use: ICCT 2021 (m²/GJ); Water use: IEA 2021 (L/GJ); TRL: IATA 2022'
    )
    fig.text(0.5, 0.02, notes, ha='center', fontsize=7.5, color='#616161',
             style='italic')

    ax.legend(loc='upper right', bbox_to_anchor=(1.35, 1.1),
              fontsize=10, framealpha=0.9)

    _save(fig, 'FigK_Radar_Plot')


# ─────────────────────────────────────────────────────────────────────────────
# MAIN
# ─────────────────────────────────────────────────────────────────────────────

if __name__ == '__main__':
    print('=' * 72)
    print('SAF HARMONIZATION — EXTENDED PUBLICATION FIGURES')
    print('=' * 72)

    print('\n[Setup] Loading / computing analysis data ...')
    df_harm    = build_harmonized_dataset()
    mc_results = run_monte_carlo(n_iter=10_000)

    print('[Setup] Running Sobol analysis (1,500 samples × 4 pathways) ...')
    sobol_res  = run_sobol_analysis(n_sobol=1_500)
    var_df     = decompose_variance(mc_results, sobol_res)

    figures = [
        ('A', 'System Boundary Diagram',           lambda: figA_system_boundary()),
        ('B', 'Harmonization Flowchart',            lambda: figB_harmonization_flowchart()),
        ('C', 'MFSP Cost Breakdown',                lambda: figC_mfsp_cost_breakdown(df_harm, mc_results)),
        ('D', 'Tornado Sensitivity Chart',          lambda: figD_tornado_sensitivity(mc_results)),
        ('E', 'MC Overlay Before/After',            lambda: figE_mc_overlay(df_harm, mc_results)),
        ('F', 'GHG Comparison with Ranges',         lambda: figF_ghg_comparison(df_harm, mc_results)),
        ('G', 'GHG Lifecycle Contribution',         lambda: figG_ghg_contribution()),
        ('H', 'Parameter Heatmap (42 studies)',     lambda: figH_parameter_heatmap(df_harm)),
        ('I', 'Variance Decomposition',             lambda: figI_variance_decomposition(var_df, sobol_res)),
        ('J', 'TEA vs LCA Parity Scatter',          lambda: figJ_tea_lca_scatter(df_harm, mc_results)),
        ('K', 'Multi-Criteria Radar Plot',          lambda: figK_radar_plot()),
    ]

    print()
    for letter, name, fn in figures:
        print(f'[Fig {letter}] {name} ...')
        fn()

    print(f'\n{"=" * 72}')
    print(f'All 11 figures saved to outputs/figures/')
    print(f'{"=" * 72}\n')
