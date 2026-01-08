#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
RIC Framework - Residual Analysis (Data minus Model)
----------------------------------------------------
Author: Gestan Morgan
Date:   January 2026
Paper:  Evidence for Saturation Limits in Dark Matter Halos (Preprint)

Description:
    This script visualizes the spectral residuals. It subtracts the 
    background power-law model from the observed flux to isolate the 
    excess signal.
    
    A distinct 'bump' in the residuals between 20-60 GeV provides 
    visual evidence of the RIC saturation prediction, distinct from 
    background noise.

Methodology:
    1. Load and stack data from 9 dSphs (Carina excluded).
    2. Bin data and calculate flux density.
    3. Fit background model to sidebands (<18 GeV, >80 GeV).
    4. Calculate Residuals = (Data - Model) / Model.
    5. Plot residuals with error bars to highlight the signal excess.
"""

import astropy.io.fits as pyfits
import matplotlib.pyplot as plt
import numpy as np
import glob
from scipy.optimize import curve_fit

# --- CONFIGURATION ---
RIC_RADIUS = 0.5           # Region of Interest radius in degrees
INPUT_SUFFIX = "_1.0.fits"
SIDEBAND_LOW = 18          # Upper limit for low-energy sideband
SIDEBAND_HIGH = 80         # Lower limit for high-energy sideband

# TARGET LIST (N=9)
# Consistent with Paper 3: Carina excluded due to contamination.
TARGETS = {
    'Bootes_1':       (210.02, 14.50),
    # 'Carina':       (100.40, -50.96), # Excluded (see Paper 3, Sec 2.2)
    'Coma_Berenices': (186.74, 23.90),
    'Draco':          (260.05, 57.91),
    'Fornax':         (39.99, -34.44),
    'Leo_1':          (152.11, 12.30),
    'Reticulum_2':    (53.92, -54.05),
    'Sculptor':       (15.03, -33.70),
    'Segue_1':        (151.76, 16.08),
    'Willman_1':      (162.34, 51.05)
}

def angular_separation(ra1, dec1, ra2, dec2):
    """ Calculates angular separation between celestial points (vectorized). """
    ra1, dec1, ra2, dec2 = map(np.radians, [ra1, dec1, ra2, dec2])
    dlon = ra2 - ra1
    dlat = dec2 - dec1
    a = np.sin(dlat/2)**2 + np.cos(dec1) * np.cos(dec2) * np.sin(dlon/2)**2
    c = 2 * np.arcsin(np.sqrt(a))
    return np.degrees(c)

def power_law(x, a, k):
    return a * x**-k

def run_residual_analysis():
    print(f"--- STARTING RIC RESIDUAL ANALYSIS (Radius: {RIC_RADIUS} deg) ---")

    all_energies = []
    processed_targets = []

    # 1. DATA LOADING
    for name, coords in TARGETS.items():
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
                target_ra, target_dec = coords
                dists = angular_separation(ph_ra, ph_dec, target_ra, target_dec)

                mask = dists < RIC_RADIUS
                valid_energies = ph_energy[mask] / 1000.0 # MeV -> GeV

                all_energies.extend(valid_energies)
                processed_targets.append(name)
                print(f"   -> {name}: {len(valid_energies)} photons loaded.")

        except Exception as e:
            print(f"❌ Error processing {name}: {e}")

    if not all_energies:
        print("Error: No data loaded.")
        return

    data = np.array(all_energies)
    print(f"\n📊 TOTAL STACK: {len(data)} photons from {len(processed_targets)} dSphs.")

    # 2. BINNING & FLUX
    bins = np.logspace(np.log10(10), np.log10(300), 14)
    counts, edges = np.histogram(data, bins=bins)
    bin_centers = (edges[:-1] * edges[1:])**0.5
    bin_width = edges[1:] - edges[:-1]

    flux = counts / bin_width
    # Poisson error propagation for flux
    flux_errors = np.sqrt(counts) / bin_width

    # 3. BACKGROUND FIT (Sidebands Only)
    mask_fit = (bin_centers < SIDEBAND_LOW) | (bin_centers > SIDEBAND_HIGH)

    try:
        popt, pcov = curve_fit(power_law, bin_centers[mask_fit], flux[mask_fit], 
                               p0=[1000, 2.0], maxfev=10000)
    except Exception as e:
        print(f"Fit Failed: {e}")
        return

    # 4. CALCULATE RESIDUALS
    model_flux = power_law(bin_centers, *popt)
    
    # Fractional Residuals: (Data - Model) / Model
    residuals = (flux - model_flux) / model_flux
    
    # Error propagation for residuals (approximate)
    res_errors = flux_errors / model_flux

    # 5. VISUALIZATION
    plt.figure(figsize=(10, 6), facecolor='white')
    
    # Reference Line (No Excess)
    plt.axhline(0, color='black', alpha=0.8, linewidth=1)

    # RIC Prediction Line
    plt.axvline(20.7, color='red', linestyle='--', alpha=0.8, linewidth=2, 
                label='RIC Prediction (20.7 GeV)')

    # Residual Data Points
    plt.errorbar(bin_centers, residuals, yerr=res_errors, fmt='o',
                 color='blue', ecolor='gray', capsize=4, markersize=8,
                 label=f'Stacked Residuals (N={len(data)})')

    # Styling
    plt.title(f'RIC Spectral Residuals (9 dSphs, Carina Excluded)', fontsize=14)
    plt.xlabel('Energy (GeV)', fontsize=12)
    plt.ylabel('Fractional Residual (Data - Model) / Model', fontsize=12)
    plt.xscale('log')
    plt.grid(True, which="both", alpha=0.2)
    plt.legend(fontsize=10)
    
    # Highlight the excess region visually
    plt.fill_between(bin_centers, 0, residuals, 
                     where=((bin_centers > 20) & (bin_centers < 60)),
                     color='red', alpha=0.1)

    outfile = '03_RIC_Residuals_Plot.png'
    plt.savefig(outfile, dpi=150)
    print(f"✅ Residual plot saved: {outfile}")

if __name__ == "__main__":
    run_residual_analysis()