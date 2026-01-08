#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
RIC Framework - Mass Invariance Check (J-Factor Regression)
-----------------------------------------------------------
Author: Gestan Morgan
Date:   January 2026
Paper:  Evidence for Saturation Limits in Dark Matter Halos (Preprint)

Description:
    This script tests the correlation between the observed gamma-ray excess
    and the gravitational mass (J-Factor) of the target dwarf galaxies.
    
    Standard WIMP Dark Matter models predict a strong positive correlation (R ~ 1.0).
    RIC predicts Mass Invariance (Saturation), resulting in no correlation (R ~ 0).

Methodology:
    1. Load J-Factors from literature (Ackermann et al. 2015, Geringer-Sameth et al. 2015).
    2. Count photons in the RIC signal window (20-60 GeV) for each target (0.5 deg ROI).
    3. Perform linear regression (SciPy) to determine Pearson R coefficient.
"""

import astropy.io.fits as pyfits
import matplotlib.pyplot as plt
import numpy as np
import glob
from scipy.stats import linregress

# --- CONFIGURATION ---
RIC_RADIUS = 0.5           # Region of Interest radius in degrees
INPUT_SUFFIX = "_1.0.fits"
SIGNAL_WINDOW = (20, 60)   # RIC Excess Window (GeV) consistent with spectral analysis

# J-FACTORS (Log10 GeV^2 cm^-5)
# Source: Ackermann et al. (2015) / Geringer-Sameth et al. (2015)
# Represents the integrated dark matter density squared along the line of sight.
TARGETS = {
    'Bootes_1':       {'coords': (210.02, 14.50),  'j_factor': 18.8},
    'Carina':         {'coords': (100.40, -50.96), 'j_factor': 18.1}, # Included for completeness, noted outlier
    'Coma_Berenices': {'coords': (186.74, 23.90),  'j_factor': 19.0},
    'Draco':          {'coords': (260.05, 57.91),  'j_factor': 18.8},
    'Fornax':         {'coords': (39.99, -34.44),  'j_factor': 18.2},
    'Leo_1':          {'coords': (152.11, 12.30),  'j_factor': 17.7},
    'Reticulum_2':    {'coords': (53.92, -54.05),  'j_factor': 18.9},
    'Sculptor':       {'coords': (15.03, -33.70),  'j_factor': 18.6},
    'Segue_1':        {'coords': (151.76, 16.08),  'j_factor': 19.5},
    'Willman_1':      {'coords': (162.34, 51.05),  'j_factor': 19.1}
}

# Optional: To reproduce the stricter subset without Carina, uncomment below:
# TARGETS.pop('Carina', None)

def angular_separation(ra1, dec1, ra2, dec2):
    """ Calculates angular separation between two celestial points. """
    ra1, dec1, ra2, dec2 = map(np.radians, [ra1, dec1, ra2, dec2])
    dlon = ra2 - ra1
    dlat = dec2 - dec1
    a = np.sin(dlat/2)**2 + np.cos(dec1) * np.cos(dec2) * np.sin(dlon/2)**2
    c = 2 * np.arcsin(np.sqrt(a))
    return np.degrees(c)

def run_mass_scaling():
    print(f"--- STARTING MASS SCALING ANALYSIS (Radius: {RIC_RADIUS} deg) ---")
    print(f"Signal Window: {SIGNAL_WINDOW[0]} - {SIGNAL_WINDOW[1]} GeV")

    j_factors = []
    photon_counts = []
    labels = []

    # 1. DATA COLLECTION
    for name, data_dict in TARGETS.items():
        filename = f"{name}{INPUT_SUFFIX}"

        if not glob.glob(filename):
            print(f"⚠️ Warning: File {filename} not found. Skipping.")
            continue

        try:
            with pyfits.open(filename) as hdul:
                data = hdul[1].data
                ph_ra = data['RA']
                ph_dec = data['DEC']
                ph_energy = data['ENERGY']

                # Spatial Filtering
                target_ra, target_dec = data_dict['coords']
                dists = angular_separation(ph_ra, ph_dec, target_ra, target_dec)
                radius_mask = dists < RIC_RADIUS

                # Energy Filtering (RIC Signal Window)
                energies_gev = ph_energy / 1000.0
                energy_mask = (energies_gev >= SIGNAL_WINDOW[0]) & (energies_gev <= SIGNAL_WINDOW[1])

                # Combined Mask
                final_mask = radius_mask & energy_mask
                count = np.sum(final_mask)

                print(f"   -> {name}: {count} photons (J={data_dict['j_factor']})")

                j_factors.append(data_dict['j_factor'])
                photon_counts.append(count)
                labels.append(name)

        except Exception as e:
            print(f"❌ Error processing {name}: {e}")

    if not j_factors:
        print("Error: No data points available.")
        return

    # 2. STATISTICAL REGRESSION
    slope, intercept, r_value, p_value, std_err = linregress(j_factors, photon_counts)

    print("\n📊 STATISTICAL RESULTS (MASS INVARIANCE CHECK)")
    print("========================================")
    print(f"Correlation Coefficient (R): {r_value:.4f}")
    print(f"Slope (m):                   {slope:.4f}")
    print(f"P-Value:                     {p_value:.6f}")
    
    # Interpretation
    if abs(r_value) < 0.2:
        print(">> RESULT: Mass Invariance CONFIRMED (R ~ 0). Supports RIC Saturation.")
    else:
        print(">> RESULT: Correlation detected. Supports WIMP/Standard Model.")
    print("========================================\n")

    # 3. VISUALIZATION
    plt.figure(figsize=(10, 7), facecolor='white')

    # Scatter Plot
    plt.scatter(j_factors, photon_counts, color='blue', s=100, edgecolors='black', 
                label=f'Observed Excess ({SIGNAL_WINDOW[0]}-{SIGNAL_WINDOW[1]} GeV)')

    # Labels
    for i, txt in enumerate(labels):
        plt.annotate(txt, (j_factors[i], photon_counts[i]),
                     xytext=(5, 5), textcoords='offset points', fontsize=9)

    # Trend Line
    x_range = np.linspace(min(j_factors)-0.2, max(j_factors)+0.2, 100)
    plt.plot(x_range, slope*x_range + intercept, 'r--', alpha=0.6, linewidth=2, 
             label=f'Regression Fit (R={r_value:.2f})')

    # Reference Line for WIMP (Hypothetical strong scaling, just for comparison if needed)
    # plt.plot(x_range, x_range * 10 - 150, 'g:', alpha=0.3, label='Expected WIMP Scaling (Slope ~ 1)')

    # Layout
    plt.title(f'RIC Mass Invariance Check (N={len(j_factors)})', fontsize=14)
    plt.xlabel('Dark Matter Density (Log10 J-Factor)', fontsize=12)
    plt.ylabel(f'Photon Counts ({SIGNAL_WINDOW[0]}-{SIGNAL_WINDOW[1]} GeV)', fontsize=12)
    plt.grid(True, alpha=0.3)
    plt.legend(fontsize=10)

    outfile = '02_RIC_Mass_Invariance.png'
    plt.savefig(outfile, dpi=150)
    print(f"✅ Plot saved to: {outfile}")

if __name__ == "__main__":
    run_mass_scaling()