# SAF Harmonization Framework - Computational Analysis

Quantifying methodological variability in sustainable aviation fuel assessments using verified data from peer-reviewed sources.

## Overview

This computational framework demonstrates that methodological choices in sustainable aviation fuel (SAF) assessments create variations of 133% in production costs and 119% in greenhouse gas emissions for identical physical systems. A three-tier harmonization framework reduces these artificial variations by 73% and 64% respectively.

## Data Sources

All calculations use publicly verifiable data from:

- **Tao et al. (2017)** Green Chemistry, DOI: [10.1039/C6GC02800D](https://doi.org/10.1039/C6GC02800D)
  - Techno-economic analysis parameters
  - Capital costs and conversion yields
  
- **Han et al. (2017)** Biotechnology for Biofuels, DOI: [10.1186/s13068-017-0698-z](https://doi.org/10.1186/s13068-017-0698-z)
  - Life cycle assessment emissions
  - Allocation methodology validation

- **GREET Model 2023** Argonne National Laboratory
  - Background emission factors
  - Available at: [greet.es.anl.gov](https://greet.es.anl.gov)

## Installation

```bash
pip install -r requirements.txt
```

## Usage

```bash
python saf_harmonization_analysis.py
```

## Outputs

The analysis generates:

1. **Excel Workbook** (`outputs/SAF_Harmonization_Analysis.xlsx`)
   - Data sources sheet with citations
   - Complete scenario results
   - Summary statistics

2. **Publication-Quality Figures** (PNG 600 DPI + PDF vector)
   - Figure 1: Main results (MFSP and GHG by scenario)
   - Figure 2: Harmonization impact (variance reduction)

## Key Findings

- **Cost variation**: $1.33 to $3.09 per gallon gasoline equivalent (133% range)
- **GHG variation**: 10.6 to 29.0 gCO₂e/MJ (119% range)
- **Harmonization effectiveness**: 73% cost variance reduction, 64% GHG variance reduction

## Methodology

The framework analyzes 14 scenarios varying:
- System boundaries (gate-to-gate, well-to-gate, well-to-pump)
- Allocation methods (energy, mass, economic, system expansion)
- Economic parameters (discount rate, feedstock cost, capacity factor)
- Environmental factors (ILUC, electricity grid intensity)

## Citation

### BibTeX
```bibtex
@software{amponsah_saf_harmonization_2026,
  author = {Amponsah, J.},
  title = {SAF Harmonization Framework - Computational Analysis},
  year = {2026},
  url = {https://github.com/joeampon/saf-harmonization},
}
```

### APA
Amponsah, J. (2026). SAF harmonization framework - computational analysis. GitHub. https://github.com/joeampon/saf-harmonization

### Chicago
Amponsah, J. "SAF Harmonization Framework - Computational Analysis." GitHub, 2026. https://github.com/joeampon/saf-harmonization.

### DOI (when available)
Request a DOI through Zenodo for formal publication citations.

**Reference the associated peer-reviewed paper when published for academic citations.**

## License

MIT License - See LICENSE file for details

## Contact & Support

For questions, bug reports, or suggestions, please open an issue on the GitHub repository. For detailed discussions, use GitHub Discussions.

## Reproducibility

All calculations are fully traceable to peer-reviewed sources. The code includes no fabricated data and all parameters can be independently verified against the cited publications.
