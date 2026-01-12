#!/usr/bin/env python3
"""
================================================================================
RIC STATISTICAL SIGNIFICANCE CALCULATOR
================================================================================

Purpose:
    Calculates the statistical significance of the gamma-ray excess in the 
    RIC signal window (18-50 GeV) using stacked Fermi-LAT data from dwarf 
    spheroidal galaxies.

Method:
    1. Load all photon energies from processed FITS files
    2. Fit a power-law background model to sidebands (<18 GeV and >60 GeV)
    3. Calculate observed vs. expected counts in the signal window
    4. Compute significance using Poisson statistics

Reference:
    Li, T.-P. & Ma, Y.-Q. (1983), ApJ 272, 317
    "Analysis methods for results in gamma-ray astronomy"

Author: Gestan Morgan
Project: RIC - Recursive Instability Collapse
Date: January 2026
================================================================================
"""

import astropy.io.fits as pyfits
import numpy as np
import glob
import os
from scipy.optimize import curve_fit
from scipy.stats import poisson


# =============================================================================
# CONFIGURATION
# =============================================================================

# Data paths for both analysis batches
PATH_BATCH_1 = "../FINAL_DATA"
PATH_BATCH_2 = "../results_2"

# Energy range for analysis (GeV)
ENERGY_MIN = 10
ENERGY_MAX = 300

# Signal window (GeV) - RIC prediction zone
SIGNAL_MIN = 18
SIGNAL_MAX = 50

# Sideband regions for background fitting (GeV)
SIDEBAND_LOW_MAX = 18
SIDEBAND_HIGH_MIN = 60

# Number of histogram bins (logarithmic spacing)
N_BINS = 15

# Significance thresholds
EVIDENCE_THRESHOLD = 3.0   # sigma
DISCOVERY_THRESHOLD = 5.0  # sigma


# =============================================================================
# ANALYSIS FUNCTIONS
# =============================================================================

def power_law(x, amplitude, index):
    """
    Power-law spectral model for background fitting.
    
    Args:
        x (array): Energy values
        amplitude (float): Normalization constant
        index (float): Spectral index
        
    Returns:
        array: Flux values
    """
    return amplitude * x**(-index)


def get_all_energies():
    """
    Loads photon energies from all processed FITS files in both batches.
    
    Returns:
        array: All photon energies in GeV
    """
    all_energies = []
    total_files = 0
    
    for path in [PATH_BATCH_1, PATH_BATCH_2]:
        files = glob.glob(os.path.join(path, '*_final.fits'))
        total_files += len(files)
        
        for filepath in files:
            with pyfits.open(filepath) as hdul:
                # Convert MeV to GeV
                energies = hdul[1].data['ENERGY'] / 1000.0
                all_energies.extend(energies)
    
    print(f"[OK] Loaded {len(all_energies)} photons from {total_files} files")
    return np.array(all_energies)


def calculate_significance():
    """
    Main analysis function. Calculates statistical significance of the 
    gamma-ray excess in the RIC signal window.
    """
    print("\n" + "="*70)
    print("RIC STATISTICAL SIGNIFICANCE ANALYSIS")
    print("="*70)
    
    # -------------------------------------------------------------------------
    # Step 1: Load data
    # -------------------------------------------------------------------------
    print("\n[1/4] Loading photon data...")
    data = get_all_energies()
    
    if len(data) == 0:
        print("[ERROR] No data found!")
        return
    
    # -------------------------------------------------------------------------
    # Step 2: Create energy histogram
    # -------------------------------------------------------------------------
    print("[2/4] Building energy spectrum...")
    
    bins = np.logspace(np.log10(ENERGY_MIN), np.log10(ENERGY_MAX), N_BINS)
    counts, edges = np.histogram(data, bins=bins)
    
    # Geometric bin centers (appropriate for log-spaced bins)
    bin_centers = np.sqrt(edges[:-1] * edges[1:])
    bin_widths = edges[1:] - edges[:-1]
    
    # Convert to flux (counts per GeV)
    flux = counts / bin_widths
    
    # -------------------------------------------------------------------------
    # Step 3: Fit background model to sidebands
    # -------------------------------------------------------------------------
    print("[3/4] Fitting background model (sidebands only)...")
    
    # Select sideband regions (exclude signal window)
    mask_sideband = (bin_centers < SIDEBAND_LOW_MAX) | (bin_centers > SIDEBAND_HIGH_MIN)
    
    # Fit power-law to sidebands
    popt, pcov = curve_fit(
        power_law, 
        bin_centers[mask_sideband], 
        flux[mask_sideband], 
        p0=[1000, 2.0]
    )
    
    print(f"      Background model: F(E) = {popt[0]:.2f} * E^(-{popt[1]:.2f})")
    
    # -------------------------------------------------------------------------
    # Step 4: Calculate significance in signal window
    # -------------------------------------------------------------------------
    print("[4/4] Calculating significance in signal window...")
    
    # Select signal region
    mask_signal = (bin_centers >= SIGNAL_MIN) & (bin_centers <= SIGNAL_MAX)
    
    # Observed counts
    observed = np.sum(counts[mask_signal])
    
    # Expected counts from background model
    expected_flux = power_law(bin_centers[mask_signal], *popt)
    expected = np.sum(expected_flux * bin_widths[mask_signal])
    
    # Excess
    excess = observed - expected
    
    # Significance (Gaussian approximation for Poisson)
    sigma = excess / np.sqrt(expected)
    
    # P-value (probability of observing this excess by chance)
    p_value = 1.0 - poisson.cdf(int(observed), expected)
    
    # -------------------------------------------------------------------------
    # Results
    # -------------------------------------------------------------------------
    print("\n" + "="*70)
    print(f"RESULTS: Signal Window {SIGNAL_MIN}-{SIGNAL_MAX} GeV")
    print("="*70)
    print(f"  Observed photons:      {int(observed)}")
    print(f"  Expected background:   {expected:.2f}")
    print(f"  Excess:                {excess:.2f} photons")
    print("-"*70)
    print(f"  SIGNIFICANCE:          {sigma:.2f} sigma")
    print(f"  P-VALUE:               {p_value:.2e}")
    print("="*70)
    
    # Interpretation
    if sigma >= DISCOVERY_THRESHOLD:
        print(f"\n[RESULT] DISCOVERY THRESHOLD EXCEEDED (>{DISCOVERY_THRESHOLD} sigma)")
        print("         This excess is statistically significant.")
    elif sigma >= EVIDENCE_THRESHOLD:
        print(f"\n[RESULT] EVIDENCE THRESHOLD EXCEEDED (>{EVIDENCE_THRESHOLD} sigma)")
        print("         Strong evidence, but below discovery threshold.")
    else:
        print(f"\n[RESULT] Below evidence threshold (<{EVIDENCE_THRESHOLD} sigma)")
        print("         Not statistically significant.")
    
    print("")
    
    return {
        'observed': int(observed),
        'expected': expected,
        'excess': excess,
        'sigma': sigma,
        'p_value': p_value
    }


# =============================================================================
# MAIN
# =============================================================================

if __name__ == "__main__":
    calculate_significance()