# SAF Harmonization Meta-Analysis Framework v2.0

**Author:** Joseph Amponsah, Iowa State University

**Email:** joeampon@iastate.edu

**Paper:** Methodology over technology: harmonizing techno-economic and life cycle assessment of sustainable aviation fuel-Biomass and Bioenergy (2026)

---

## Overview

This repository contains the full computational framework for the systematic harmonization of 48 peer-reviewed SAF (sustainable aviation fuel) TEA and LCA studies across four production pathways:

| Pathway | n  | Description                                          |
| ------- | -- | ---------------------------------------------------- |
| ATJ     | 8  | Alcohol-to-Jet (corn stover, sugarcane, wheat straw) |
| HEFA    | 10 | Hydroprocessed Esters and Fatty Acids                |
| FT-SPK  | 14 | Fischer-Tropsch Synthetic Paraffinic Kerosene        |
| PtL     | 16 | Power-to-Liquid (electrolysis + DAC + FT)            |

---

## Harmonization Protocol

Five-step protocol applied uniformly to all 48 studies:

| Step | Action                                                  | Reference                   |
| ---- | ------------------------------------------------------- | --------------------------- |
| 1    | Convert MFSP → 2023 USD/GGE (CPI + unit conversion)    | US BLS CPI-U                |
| 2    | Re-weight allocation → energy basis                    | ISO 14044 §4.3.4.2         |
| 3    | Adjust system boundary → Well-to-Wake (+3.0 gCO₂e/MJ) | ICAO CORSIA Doc 10164       |
| 4    | Remove ILUC if study included it                        | ICAO CORSIA ILUC table      |
| 5    | CRF normalization (10% DR / 30-yr / 90% CF)             | NREL design case convention |

**Key fix vs v1.0:** Step 5 now uses pathway-specific CAPEX fractions (ATJ=0.41, HEFA=0.13, FT-SPK=0.55, PtL=0.19) instead of a uniform 0.40. This prevents over-correction in feedstock-dominated pathways like HEFA.

---

## Project Structure

```
saf_harmonization/
├── main.py                          ← Entry point
├── config.py                        ← Constants, escalation tables, paths
├── requirements.txt
├── data/
│   ├── literature_database.py       ← 48 studies with full metadata + DOIs
│   ├── parameter_distributions.py   ← Monte Carlo distributions + sources
│   └── generate_excel.py            ← 9-sheet Excel workbook generator
├── harmonization/
│   └── engine.py                    ← Five-step harmonization protocol
├── models/
│   └── pathway_models.py            ← ATJ, HEFA, FT-SPK, PtL TEA+LCA models
├── analysis/
│   ├── monte_carlo.py               ← Monte Carlo simulation (10,000 iter)
│   ├── sobol_analysis.py            ← Jansen (1999) Sobol estimator
│   └── variance_decomposition.py    ← Variance decomp + external validation
├── visualization/
│   └── figures.py                   ← 16 publication figures (PNG + PDF)
└── outputs/
    ├── figures/                     ← Generated figures
    └── SAF_MetaAnalysis_Harmonization.xlsx
```

---

## Quick Start

```bash
# 1. Clone repository
git clone https://github.com/joeampon/saf-harmonization-main.git
cd saf-harmonization

# 2. Install dependencies
pip install -r requirements.txt

# 3. Quick test run (~2 min)
python main.py --fast

# 4. Full analysis (~25 min)
python main.py
```

---

## Outputs

### Figures (outputs/figures/)

| Figure | Content                                            |
| ------ | -------------------------------------------------- |
| Fig1   | Raw vs harmonized literature scatter (MFSP vs GHG) |
| Fig2   | Harmonization impact before/after per pathway      |
| Fig3   | Monte Carlo violin distributions                   |
| Fig4   | Variance decomposition (3 panels)                  |
| FigA   | Well-to-Wake system boundary schematic             |
| FigB   | Five-step harmonization flowchart                  |
| FigC   | MFSP cost component breakdown                      |
| FigD   | Tornado OAT sensitivity chart                      |
| FigE   | MC KDE overlay vs literature                       |
| FigF   | GHG distributions vs regulatory thresholds         |
| FigG   | GHG component breakdown                            |
| FigH   | Sobol S1 heatmap                                   |
| FigI   | Extended variance decomposition (6 panels)         |
| FigJ   | TEA vs LCA trade-off scatter                       |
| FigK   | Multi-dimensional radar plot                       |
| FigL   | Pareto Sobol sensitivity chart                     |

All figures saved at 300 DPI (PNG) and as vector PDF with embedded fonts (Elsevier requirement).

### Excel Workbook (9 sheets)

1. Literature Database — studies with DOIs
2. Harmonized Values — step-by-step per study
3. Parameter Distributions — Monte Carlo ranges with sources
4. Monte Carlo Summary — P5/P25/median/mean/P75/P95/Std/CV
5. Sobol Indices — S1 and ST per parameter
6. Variance Decomposition — methodological vs technical split
7. External Validation — 10 independent 2025 studies
8. Key Findings — headline numbers
9. References — complete reference list

---

## Key Results

- Methodology creates **133% cost variation** ($1.33–$3.09/GGE) for identical corn stover ATJ systems
- Methodology creates **174% GHG variation** (10.6–29.0 gCO₂e/MJ) for identical systems
- Three-tier framework reduces cost CV by **63%** and GHG CV by **85%**
- All 48 harmonized studies show >50% GHG reduction vs petroleum jet (89 gCO₂e/MJ)
- Only 3 HEFA studies (palm/soy + ILUC) exceed the EU RED III threshold (31.15 gCO₂e/MJ)

---

## Reference Basis

| Parameter         | Value              | Source                       |
| ----------------- | ------------------ | ---------------------------- |
| Petroleum jet GHG | 89 gCO₂e/MJ       | ICAO CORSIA Doc 10164 (2022) |
| Cost year         | 2023 USD           | CEPCI 2023 = 798             |
| Discount rate     | 10% real           | NREL Nth-plant convention    |
| Plant lifetime    | 30 years           | NREL convention              |
| Capacity factor   | 90%                | Industry standard            |
| System boundary   | Well-to-Wake       | ICAO CORSIA                  |
| Allocation        | Energy (ISO 14044) | ISO 14044 §4.3.4.2          |

---

## Citation

```bibtex
@article{amponsah2026saf,
  title   = {Harmonized techno-economic and life cycle assessment of
             corn stover alcohol-to-sustainable aviation fuel},
  author  = {Amponsah, Joseph and Ansah, Kingsford Kweku and
             Asare, Daniel Baffour and Siameh, Fredrick Ofori and
             Owusu, Emmanuel and Essumang, Albert Adjekum and Morgan, Louis},
  journal = {Biomass and Bioenergy},
  year    = {2026},
  doi     = {10.1016/j.biombioe.2026.109261}
}
```

---

## License

MIT License — see LICENSE file.
