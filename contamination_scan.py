"""
RIC Project: Fermi-LAT 4FGL Point Source Contamination Scanner
Author: Gestan Morgan (RIC Framework)
Date: January 12, 2026

Description:
This script performs a systematic cross-match between target Dwarf Spheroidal 
Galaxies (dSphs) and the Fermi-LAT 14-year Source Catalog (4FGL-DR4). 
It identifies potential gamma-ray point source contamination within a 
defined Region of Interest (ROI) to ensure the integrity of Dark Matter 
signal extraction.
"""

import numpy as np
from astropy.io import fits
from astropy.coordinates import SkyCoord
import astropy.units as u
import os
import urllib.request
import json

# Target Dictionary: 25 Dwarf Spheroidal Galaxies (dSphs)
# Coordinates (J2000) and literature J-Factors
TARGETS = {
    "Leo_II":           {"ra": 168.37, "dec": 22.15,  "j_factor": 17.5},
    "Canes_Venatici_I": {"ra": 202.53, "dec": 33.56,  "j_factor": 17.2},
    "Ursa_Major_I":     {"ra": 158.00, "dec": 51.92,  "j_factor": 17.2},
    "Ursa_Major_II":    {"ra": 132.87, "dec": 63.13,  "j_factor": 17.9},
    "Canes_Venatici_II":{"ra": 194.29, "dec": 34.32,  "j_factor": 17.1},
    "Leo_IV":           {"ra": 173.24, "dec": -0.55,  "j_factor": 16.3},
    "Hercules":         {"ra": 247.76, "dec": 12.79,  "j_factor": 16.9},
    "Segue_2":          {"ra": 35.82,  "dec": 20.18,  "j_factor": 16.2},
    "Triangulum_II":    {"ra": 33.32,  "dec": 36.18,  "j_factor": 19.1},
    "Hydra_II":         {"ra": 185.43, "dec": -31.98, "j_factor": 16.2},
    "Eridanus_II":      {"ra": 56.09,  "dec": -43.54, "j_factor": 15.6},
    "Aquarius_III":     {"ra": 331.06, "dec": -10.58, "j_factor": 16.0},
    "Ursa_Minor":       {"ra": 227.28, "dec": 67.22,  "j_factor": 18.8},
    "Segue_1":          {"ra": 151.76, "dec": 16.08,  "j_factor": 19.5},
    "Willman_1":        {"ra": 162.34, "dec": 51.05,  "j_factor": 19.1},
    "Draco":            {"ra": 260.05, "dec": 57.91,  "j_factor": 18.8},
    "Coma_Berenices":   {"ra": 186.74, "dec": 23.90,  "j_factor": 19.0},
    "Reticulum_2":      {"ra": 53.92,  "dec": -54.05, "j_factor": 18.9},
    "Bootes_1":         {"ra": 210.02, "dec": 14.50,  "j_factor": 18.8},
    "Sculptor":         {"ra": 15.04,  "dec": -33.71, "j_factor": 18.6},
    "Fornax":           {"ra": 40.00,  "dec": -34.45, "j_factor": 18.2},
    "Carina":           {"ra": 100.40, "dec": -50.97, "j_factor": 18.1},
    "Sextans":          {"ra": 153.26, "dec": -1.02,  "j_factor": 18.4},
    "Leo_I":            {"ra": 152.12, "dec": 12.31,  "j_factor": 17.7},
    "Tucana_2":         {"ra": 342.98, "dec": -58.57, "j_factor": 19.3}
}

# Configuration
SEARCH_RADIUS = 0.5  # degrees (ROI size)
CATALOG_URL = "https://fermi.gsfc.nasa.gov/ssc/data/access/lat/14yr_catalog/gll_psc_v32.fit"
CATALOG_FILE = "gll_psc_v32.fit"

def download_catalog():
    """Downloads the Fermi 4FGL-DR4 catalog if not present locally."""
    if os.path.exists(CATALOG_FILE):
        print(f"[INFO] Catalog file found: {CATALOG_FILE}")
        return True
    
    print("[INFO] Downloading 4FGL catalog from NASA Fermi servers...")
    try:
        urllib.request.urlretrieve(CATALOG_URL, CATALOG_FILE)
        print("[SUCCESS] Download complete.")
        return True
    except Exception as e:
        print(f"[ERROR] Download failed: {e}")
        return False

def load_catalog():
    """Loads FITS catalog and extracts astrophysical parameters."""
    print(f"[INFO] Opening catalog: {CATALOG_FILE}")
    with fits.open(CATALOG_FILE) as hdul:
        data = hdul[1].data
        catalog = {
            'name': data['Source_Name'],
            'ra': data['RAJ2000'],
            'dec': data['DEJ2000'],
            'class': data['CLASS1'],
            'flux': data['Flux1000'],
            'significance': data['Signif_Avg']
        }
        print(f"[INFO] Successfully loaded {len(catalog['name'])} sources.")
        return catalog

def interpret_class(class_code):
    """Maps 4FGL classification codes to descriptive strings."""
    class_dict = {
        'bll': 'BL Lacertae object', 'fsrq': 'Flat-spectrum radio quasar',
        'bcu': 'Blazar candidate of uncertain type', 'psr': 'Pulsar',
        'pwn': 'Pulsar wind nebula', 'snr': 'Supernova remnant',
        'agn': 'Other non-blazar active galactic nucleus', 
        'rdg': 'Radio galaxy', 'sbg': 'Starburst galaxy',
        'sey': 'Seyfert galaxy', 'glc': 'Globular cluster',
        'spp': 'Unidentified (potential PSR)', 'unk': 'Unknown'
    }
    code = str(class_code).strip().lower()
    return class_dict.get(code, f'Unclassified ({class_code})')

def evaluate_contamination(contamination_list):
    """Analyzes identified sources and provides inclusion recommendations."""
    critical_types = ['bll', 'fsrq', 'bcu', 'psr', 'agn']
    
    for source in contamination_list:
        code = str(source['class']).strip().lower()
        if source['separation'] < 0.2:
            return "EXCLUDE: Source proximity to target center (< 0.2 deg)"
        if code in critical_types and source['flux'] > 1e-10:
            return "EXCLUDE: Significant gamma-ray emitter detected in ROI"
        if source['significance'] > 10.0:
            return "EXCLUDE: High-significance source detected in ROI"
            
    return "CAUTION: Low-level contamination; manual verification recommended"

def perform_scan(catalog):
    """Executes spatial cross-matching between targets and 4FGL sources."""
    cat_coords = SkyCoord(ra=catalog['ra']*u.deg, dec=catalog['dec']*u.deg, frame='icrs')
    results = {}

    print("-" * 75)
    print(f"{'Target Name':<20} | {'Status':<12} | {'Recommendation'}")
    print("-" * 75)

    for name, info in TARGETS.items():
        target_coord = SkyCoord(ra=info['ra']*u.deg, dec=info['dec']*u.deg, frame='icrs')
        separations = target_coord.separation(cat_coords).deg
        nearby_indices = np.where(separations < SEARCH_RADIUS)[0]

        if len(nearby_indices) == 0:
            status, rec = "CLEAN", "INCLUDE"
            results[name] = {"status": status, "contamination": [], "recommendation": rec}
        else:
            contamination_list = []
            for idx in nearby_indices:
                contamination_list.append({
                    'name': catalog['name'][idx],
                    'separation': separations[idx],
                    'class': catalog['class'][idx],
                    'flux': catalog['flux'][idx],
                    'significance': catalog['significance'][idx]
                })
            status = "CONTAMINATED"
            rec = evaluate_contamination(contamination_list)
            results[name] = {"status": status, "contamination": contamination_list, "recommendation": rec}

        print(f"{name:<20} | {status:<12} | {rec}")
    
    return results

def export_data(results):
    """Saves scan results to JSON for documentation purposes."""
    output_file = "contamination_analysis.json"
    serializable_data = {}
    for name, data in results.items():
        serializable_data[name] = {
            "status": data["status"],
            "recommendation": data["recommendation"],
            "sources": [
                {
                    "name": str(s['name']),
                    "separation_deg": float(s['separation']),
                    "type": interpret_class(s['class']),
                    "flux": float(s['flux']) if not np.isnan(s['flux']) else None,
                    "sig": float(s['significance'])
                } for s in data["contamination"]
            ]
        }
    with open(output_file, 'w') as f:
        json.dump(serializable_data, f, indent=4)
    print(f"\n[INFO] Data exported to {output_file}")

if __name__ == "__main__":
    print("=" * 75)
    print("   RIC FRAMEWORK - 4FGL SOURCE CONTAMINATION SCANNER")
    print("=" * 75)
    if download_catalog():
        cat_data = load_catalog()
        scan_results = perform_scan(cat_data)
        export_data(scan_results)
        print("-" * 75)
        print("Analysis complete.")