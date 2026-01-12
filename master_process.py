#!/usr/bin/env python3
"""
================================================================================
RIC FERMI-LAT DATA PROCESSING PIPELINE
================================================================================

Purpose:
    Processes raw Fermi-LAT photon data for dwarf spheroidal galaxy analysis.
    Applies strict geometric filtering (0.5 degree ROI) to isolate core halo 
    emission and minimize contamination from nearby sources.

Pipeline Steps:
    1. gtselect  - Spatial, temporal, and energy filtering
    2. Header fix - Ensures correct coordinate metadata in FITS files
    3. gtmktime  - Good Time Interval (GTI) selection based on data quality

Requirements:
    - Fermi Science Tools (fermitools) installed and in PATH
    - Input files: L*PH*.fits (photon data), L*SC*.fits (spacecraft data)
    - astropy

Output:
    - {target_name}_final.fits in ../results_2/ directory

Author: Gestan Morgan
Project: RIC - Recursive Instability Collapse
Date: January 2026
================================================================================
"""

import os
import glob
import sys
from astropy.io import fits


# =============================================================================
# CONFIGURATION
# =============================================================================

# Region of Interest radius (degrees) - strict filtering for core halo emission
ROI_RADIUS = 0.5

# Energy range (MeV) - 10 GeV to 300 GeV
ENERGY_MIN = 10000
ENERGY_MAX = 300000

# Time range (MET seconds) - August 2008 to January 2025
TIME_MIN = 239557417
TIME_MAX = 757555200

# Maximum zenith angle (degrees)
ZENITH_MAX = 180

# Output directory
OUTPUT_DIR = "../../results_2/"


# =============================================================================
# HELPER FUNCTIONS
# =============================================================================

def run_command(cmd):
    """
    Executes a shell command and exits on failure.
    
    Args:
        cmd (str): Command to execute
        
    Raises:
        SystemExit: If command returns non-zero exit code
    """
    print(f"[CMD] {cmd}")
    if os.system(cmd) != 0:
        print(f"[ERROR] Command failed: {cmd}")
        sys.exit(1)


def find_input_files():
    """
    Automatically locates spacecraft and photon files in current directory.
    
    Returns:
        tuple: (spacecraft_file, list_of_photon_files)
        
    Raises:
        SystemExit: If required files are not found
    """
    sc_files = glob.glob('*_SC*.fits')
    ph_files = [f for f in glob.glob('L*.fits') 
                if 'SC' not in f 
                and 'filtered' not in f 
                and 'final' not in f]
    
    if not sc_files or not ph_files:
        print("[ERROR] Required files not found!")
        print("        Need: L*PH*.fits (photon data)")
        print("        Need: L*SC*.fits (spacecraft data)")
        sys.exit(1)
    
    return sc_files[0], ph_files


def get_target_coordinates():
    """
    Prompts user for target coordinates and name.
    
    Returns:
        tuple: (ra, dec, target_name)
    """
    print("\n--- TARGET INFORMATION ---")
    ra = input("RA (degrees):  ").strip()
    dec = input("DEC (degrees): ").strip()
    name = input("Target name:   ").strip()
    return ra, dec, name


def fix_fits_header(filename, ra, dec):
    """
    Updates FITS header with correct coordinate metadata.
    
    Args:
        filename (str): Path to FITS file
        ra (str): Right Ascension in degrees
        dec (str): Declination in degrees
    """
    try:
        with fits.open(filename, mode='update') as hdul:
            header = hdul['EVENTS'].header
            header['RA_OBJ'] = float(ra)
            header['DEC_OBJ'] = float(dec)
            header['RADIUS'] = ROI_RADIUS
            hdul.flush()
        print("[OK] Header metadata updated")
    except Exception as e:
        print(f"[WARNING] Header fix skipped: {e}")


# =============================================================================
# MAIN PIPELINE
# =============================================================================

def main():
    """
    Main processing pipeline for Fermi-LAT dSph analysis.
    """
    print("\n" + "="*70)
    print("RIC FERMI-LAT PROCESSING PIPELINE")
    print(f"ROI Radius: {ROI_RADIUS} degrees | Energy: {ENERGY_MIN/1000:.0f}-{ENERGY_MAX/1000:.0f} GeV")
    print("="*70)
    
    # Locate input files
    sc_file, ph_files = find_input_files()
    print(f"\n[OK] Found spacecraft file: {sc_file}")
    print(f"[OK] Found {len(ph_files)} photon file(s)")
    
    # Get target information
    ra, dec, target_name = get_target_coordinates()
    
    # Create event list file for gtselect
    with open("events.txt", "w") as f:
        for ph in ph_files:
            f.write(ph + "\n")
    
    # -------------------------------------------------------------------------
    # Step 1: gtselect - Spatial, temporal, and energy filtering
    # -------------------------------------------------------------------------
    print(f"\n[1/3] Running gtselect (ROI: {ROI_RADIUS} deg)...")
    
    gtselect_cmd = (
        f"gtselect "
        f"infile=@events.txt "
        f"outfile={target_name}_filtered.fits "
        f"ra={ra} dec={dec} rad={ROI_RADIUS} "
        f"tmin={TIME_MIN} tmax={TIME_MAX} "
        f"emin={ENERGY_MIN} emax={ENERGY_MAX} "
        f"zmax={ZENITH_MAX}"
    )
    run_command(gtselect_cmd)
    
    # -------------------------------------------------------------------------
    # Step 2: Fix FITS header metadata
    # -------------------------------------------------------------------------
    print("\n[2/3] Fixing FITS header...")
    fix_fits_header(f"{target_name}_filtered.fits", ra, dec)
    
    # -------------------------------------------------------------------------
    # Step 3: gtmktime - Good Time Interval selection
    # -------------------------------------------------------------------------
    print("\n[3/3] Running gtmktime (GTI selection)...")
    
    gtmktime_cmd = (
        f"gtmktime "
        f"scfile={sc_file} "
        f"filter='DATA_QUAL>0 && LAT_CONFIG==1' "
        f"roicut=no "
        f"evfile={target_name}_filtered.fits "
        f"outfile={target_name}_final.fits"
    )
    run_command(gtmktime_cmd)
    
    # -------------------------------------------------------------------------
    # Step 4: Move results and cleanup
    # -------------------------------------------------------------------------
    print(f"\n[4/4] Moving results to {OUTPUT_DIR}")
    run_command(f"mv {target_name}_final.fits {OUTPUT_DIR}")
    
    # Cleanup intermediate files
    os.system(f"rm {target_name}_filtered.fits events.txt 2>/dev/null")
    
    # -------------------------------------------------------------------------
    # Done
    # -------------------------------------------------------------------------
    print("\n" + "="*70)
    print(f"[COMPLETE] {target_name} processed successfully")
    print(f"Output: {OUTPUT_DIR}{target_name}_final.fits")
    print("="*70 + "\n")


if __name__ == "__main__":
    main()