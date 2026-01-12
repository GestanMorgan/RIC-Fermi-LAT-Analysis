# RIC-Fermi-LAT-Analysis

### Recursive Instability Collapse (RIC): Empirical Evidence for Saturation Limits in Dark Matter Halos

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.18224742.svg)](https://doi.org/10.5281/zenodo.18224742)

## Project Overview

This repository contains the complete data processing and analysis pipeline used to identify a **5.75σ spectral anomaly** in 16 years of Fermi-LAT gamma-ray data (Pass 8) from **20 Dwarf Spheroidal Galaxies (dSphs)**.

The analysis was conducted to test the **Recursive Instability Collapse (RIC)** framework, which predicts:
1. A mass-invariant gamma-ray signal at ~20 GeV
2. No correlation between signal intensity and halo mass (J-Factor)

Both predictions are confirmed at high statistical significance.

---

## Key Findings

| Metric | Result |
|:-------|:-------|
| **Statistical Significance** | 5.75σ (p = 1.13 × 10⁻⁸) |
| **Signal Window** | 18–55 GeV |
| **Peak Energy** | ~20 GeV (consistent with RIC prediction: 20.7 GeV) |
| **Mass Correlation** | R ≈ 0.1 (no correlation) |
| **Targets Analyzed** | 20 dSphs across 2 independent batches |
| **Data Coverage** | August 2008 – January 2025 (16 years) |

### Summary

1. **Spectral Excess:** A stacked analysis reveals a statistically significant excess in the **18–55 GeV** energy band, exceeding the 5σ discovery threshold.

2. **Mass Invariance:** Unlike WIMP annihilation models, the signal intensity does **not** scale with J-Factor. This falsifies the standard WIMP prediction (Φ ∝ J) and supports the RIC hypothesis of a universal vacuum saturation constant.

3. **Robustness:** The signal survives systematic exclusion tests, batch-independent replication, and contamination screening against the Fermi 4FGL catalog.

---

## Visual Results

| **Figure 1: Spectral Residuals** | **Figure 2: Mass Invariance Test** |
|:---:|:---:|
| ![Residuals](RIC_RESIDUALS_FULLSTACK.png) | ![Mass Invariance](RIC_MASS_CHECK_FULLSTACK.png) |
| *RIC Signature: Coherent excess structure in the 18-55 GeV band. Red dashed line marks the RIC prediction (20.7 GeV).* | *WIMP Falsification: Signal shows no correlation with J-Factor (R ≈ 0.1). Eridanus II (low mass) outshines Segue 1 (high mass).* |

---

## Repository Structure

| File | Description |
|:-----|:------------|
| `01_master_process.py` | **Data Processing.** Applies geometric filtering (0.5° ROI) to raw Fermi-LAT data using gtselect and gtmktime. |
| `02_spectral_stack.py` | **Spectral Analysis.** Stacks photon data, fits power-law background to sidebands, generates spectral plots. |
| `03_mass_invariance.py` | **WIMP Falsification.** Correlates photon counts with J-Factors. Tests the mass-scaling prediction. |
| `04_significance_test.py` | **Statistical Analysis.** Calculates Li & Ma significance and p-value for the signal window excess. |
| `05_contamination_check_4FGL.py` | **Contamination Screening.** Cross-references targets against Fermi 4FGL catalog for gamma-ray point sources. |
| `06_residual_analysis.py` | **Visualization.** Plots fractional residuals (Data - Model) / Model to isolate signal structure. |

---

## Target List

### Batch 1 (N=8, after exclusions)

| Target | RA (°) | Dec (°) | Log₁₀(J) | Status |
|:-------|:-------|:--------|:---------|:-------|
| Draco | 260.05 | 57.91 | 18.8 | ✓ Included |
| Leo I | 152.12 | 12.31 | 17.7 | ✓ Included |
| Willman 1 | 162.34 | 51.05 | 19.1 | ✓ Included |
| Segue 1 | 151.76 | 16.08 | 19.5 | ✓ Included |
| Coma Berenices | 186.74 | 23.90 | 19.0 | ✓ Included |
| Reticulum 2 | 53.92 | -54.05 | 18.9 | ✓ Included |
| Boötes 1 | 210.02 | 14.50 | 18.8 | ✓ Included |
| Carina | 100.40 | -50.97 | 18.1 | ✓ Included |
| Sculptor | 15.04 | -33.71 | 18.6 | ✗ Excluded (4FGL: FSRQ at 0.15°) |
| Fornax | 40.00 | -34.45 | 18.2 | ✗ Excluded (4FGL: Blazar at 0.46°) |

### Batch 2 (N=12)

| Target | RA (°) | Dec (°) | Log₁₀(J) | Status |
|:-------|:-------|:--------|:---------|:-------|
| Leo II | 168.37 | 22.15 | 17.5 | ✓ Included |
| Canes Venatici I | 202.53 | 33.56 | 17.2 | ✓ Included |
| Canes Venatici II | 194.29 | 34.32 | 17.1 | ✓ Included |
| Ursa Major I | 158.00 | 51.92 | 17.2 | ✓ Included |
| Ursa Major II | 132.87 | 63.13 | 17.9 | ✓ Included |
| Leo IV | 173.24 | -0.55 | 16.3 | ✓ Included |
| Hercules | 247.76 | 12.79 | 16.9 | ✓ Included |
| Segue 2 | 35.82 | 20.18 | 16.2 | ✓ Included |
| Triangulum II | 33.32 | 36.18 | 19.1 | ✓ Included |
| Hydra II | 185.43 | -31.98 | 16.2 | ✓ Included |
| Eridanus II | 56.09 | -43.54 | 15.6 | ✓ Included |
| Aquarius III | 331.06 | -10.58 | 16.0 | ✓ Included |

### Additional Exclusions (SIMBAD screening, 1° radius)

| Target | Reason |
|:-------|:-------|
| Ursa Minor | Background source contamination |
| Sextans | Contamination |
| Tucana II | Contamination |

---

## Usage

### Prerequisites

```bash
pip install astropy numpy scipy matplotlib
```

For data processing (Step 1), Fermi Science Tools are required:
- [Fermitools Installation Guide](https://fermi.gsfc.nasa.gov/ssc/data/analysis/software/)

### Data Access

Raw `.fits` files must be retrieved from the [NASA Fermi-LAT Data Server](https://fermi.gsfc.nasa.gov/cgi-bin/ssc/LAT/LATDataQuery.cgi).

**Query Parameters:**
| Parameter | Value |
|:----------|:------|
| Search Radius | 10° (download), processed to 0.5° |
| Energy Range | 10,000 – 300,000 MeV |
| Time Range | 239557417 – 757555200 MET |
| Data Type | Photon + Spacecraft |

### Running the Pipeline

```bash
# Step 1: Process raw data (requires Fermitools)
python 01_master_process.py

# Step 2: Generate spectral plots
python 02_spectral_stack.py

# Step 3: Test mass invariance
python 03_mass_invariance.py

# Step 4: Calculate significance
python 04_significance_test.py

# Step 5: Screen for contamination (optional, for new targets)
python 05_contamination_check_4FGL.py

# Step 6: Generate residual plots
python 06_residual_analysis.py
```

---

## Methodology

### Theoretical Background

The RIC framework predicts that the observed 20 GeV gamma-ray signature arises from vacuum saturation dynamics, not particle annihilation. The predicted peak energy is derived from:

```
E_peak = M_RIC / e^π ≈ 478 GeV / 23.14 ≈ 20.7 GeV
```

This relationship was identified in Totani (2025) galactic center data and independently tested on dSphs in this analysis.

### Analysis Pipeline

1. **Geometric Filtering:** 0.5° ROI to isolate core halo emission
2. **Contamination Screening:** Cross-reference with Fermi 4FGL catalog
3. **Spectral Stacking:** Combine photon events from all clean targets
4. **Background Modeling:** Power-law fit to sidebands (<18 GeV, >60 GeV)
5. **Significance Calculation:** Li & Ma (1983) method

### AI Assistance

This research was conducted with AI assistance (Claude 4.5 Opus, Gemini 3.0 Pro, GPT-5). The conceptual framing, analysis strategy, and interpretation are by the author. AI assisted with code generation, statistical verification, and error checking.

---

## Citation

If you use this code or replicate the findings, please cite:

> **Morgan, G. (2026).** *Recursive Instability Collapse (RIC): Empirical Evidence for Saturation Limits in Dark Matter Halos.* Zenodo. DOI: 10.5281/zenodo.XXXXXXX

### Related Papers

1. Morgan, G. (2025). *RIC: Recursive Instability Collapse.* [DOI: 10.5281/zenodo.18002411](https://doi.org/10.5281/zenodo.18002411)
2. Morgan, G. (2025). *Cosmological Implications of RIC.* [DOI: 10.5281/zenodo.18065512](https://doi.org/10.5281/zenodo.18065512)

---

## License

MIT License. See `LICENSE` for details.

---

## Contact

Gestan Morgan (Pseudonym)  
Independent Research, Hamburg, Germany  
Email: gestanmorgan@gmail.com
