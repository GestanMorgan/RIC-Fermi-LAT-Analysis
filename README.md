# RIC-Fermi-LAT-Analysis

# RIC-Fermi-Analysis
### Recursive Instability Collapse (RIC): Empirical Evidence for Saturation Limits in Dark Matter Halos

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.XXXXXX.svg)](https://doi.org/10.5281/zenodo.XXXXXX)

## Project Overview
This repository contains the data processing pipeline used to identify a **3.82σ spectral anomaly** in Fermi-LAT gamma-ray data (Pass 8) from 9 Dwarf Spheroidal Galaxies (dSphs).

The analysis was conducted to test the **Recursive Instability Collapse (RIC)** framework, which predicts a mass-invariant saturation signal due to thermodynamic information limits (Bremermann Limit) in galactic cores.

### Key Findings

1.  **Spectral Excess:** A stacked analysis reveals a statistically significant excess in the **20-60 GeV** energy band ($3.82\sigma$ local significance).
2.  **Mass Invariance:** Unlike WIMP annihilation models, the signal intensity does not scale with the J-Factor (Gravitational Mass). Regression analysis yields a correlation coefficient of $R \approx -0.10$, supporting a saturation-based mechanism over particle annihilation.
3.  **Robustness:** The signal persists after excluding potential outliers (e.g., Carina) and applying strict geometric filtering ($0.5^\circ$ ROI).

## Visual Results

| **Figure 1: Mass Invariance (WIMP Falsification)** | **Figure 2: Spectral Stack (RIC Signature)** |
|:---:|:---:|
| ![Mass Scaling Plot](02_RIC_MASS_SCALING.png) | ![Spectral Stack Plot](RIC_STACK_SPECTRUM.png) |
| *Analysis showing zero correlation ($R \approx -0.10$) between Halo Mass and Photon Flux. Validated via Script 02.* | *Significant excess ($3.82\sigma$) observed in the 20-60 GeV band (Clean Sample). Validated via Script 01.* |

## Repository Structure

The analysis is modularized into three distinct steps for transparency and reproducibility:

| File | Description |
| :--- | :--- |
| `01_spectral_stack.py` | **Primary Analysis.** Stacks photon data from 9 dSphs, fits a background model (sidebands), and calculates the local significance ($3.82\sigma$) of the excess. |
| `02_mass_invariance.py` | **Falsification Test.** Correlates photon flux with J-Factors (Gravitational Mass). Result: $R \approx -0.10$, falsifying the WIMP hypothesis. |
| `03_residual_analysis.py` | **Visualization.** Plots the explicit spectral residuals `(Data - Model) / Model` to visualize the structural distinctness of the signal excess. |

*(Note: Raw .fits files must be placed in the root directory for these scripts to execute.)*

## Usage

### Prerequisites
* Python 3.8+
* `astropy`, `numpy`, `scipy`, `matplotlib`

### Data Access (Reproduction)
Due to repository size limits, raw `.fits` files must be retrieved from the [NASA Fermi-LAT Data Server](https://fermi.gsfc.nasa.gov/cgi-bin/ssc/LAT/LATDataQuery.cgi).

**Query Parameters used:**
* **Search Radius:** 10 degrees (download), analyzed at 1.0 degree.
* **Energy Range:** 1,000 - 300,000 MeV (1-300 GeV).
* **Target List:** See table below.

### Target List & Parameters
The following Dwarf Spheroidal Galaxies (dSphs) were selected based on the standard criteria for Dark Matter searches.

| Target Name | RA (deg) | Dec (deg) | Log10(J) | Status |
| :--- | :--- | :--- | :--- | :--- |
| **Draco** | 260.05 | +57.91 | 18.8 | Included |
| **Ursa Minor** | 227.28 | +67.23 | 18.8 | *Excluded (Background Source)* |
| **Carina** | 100.40 | -50.96 | 18.1 | *Excluded (Outlier)* |
| **Fornax** | 39.99 | -34.44 | 18.2 | Included |
| **Sculptor** | 15.03 | -33.70 | 18.6 | Included |
| **Leo I** | 152.11 | +12.30 | 17.8 | Included |
| **Sextans** | 153.26 | -1.61 | 17.9 | *Excluded (Contamination)* |
| **Leo II** | 168.37 | +22.15 | 17.6 | Included |
| **Canes Venatici I** | 202.01 | +33.56 | 17.4 | Included |
| **Willman 1** | 162.34 | +51.05 | 19.1 | Included |
| **Segue 1** | 151.76 | +16.08 | 19.5 | Included |

### Running the Pipeline
Place the downloaded `.fits` files in the root directory (Naming: `[GalaxyName]_1.0.fits`) and execute:

```bash
python 01_spectral_stack.py
python 02_mass_invariance.py
python 03_residual_analysis.py



## Methodology & AI-Orchestration

This research utilizes a **Recursive Emergence** methodology.
* **Theory:** The RIC framework posits that Dark Matter halos are subject to information-theoretic saturation limits ($K_{vac}$).
* **Orchestration:** The theoretical framework and analytical strategy were directed by **Gestan Morgan**.
* **Computation:** Advanced AI models served as constructive partners for code generation, error checking, and statistical formalization, ensuring a rigorous "human-in-the-loop" verification process.

## Citation
If you use this code or replicate the findings, please cite the associated preprint:

> **Morgan, G. (2026).** *Recursive Instability Collapse: The Bremermann Limit as a Cosmological Regulator.* Zenodo. DOI: [Insert Link]


