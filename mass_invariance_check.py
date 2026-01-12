import numpy as np
import glob
import matplotlib.pyplot as plt
from astropy.io import fits
from scipy.stats import pearsonr
import os

# ------------------------------------------------------------
# RIC DATABASE (BATCH 2)
# Log10 J-factors (units: GeV^2 cm^-5)
# Used here as a proxy for halo mass / information density M
#
# IMPORTANT:
# Galaxy names MUST exactly match the FITS filenames
# (everything before "_final.fits")

# This script tests the mass invariance prediction of the
#Recursive Instability Collapse (RIC) framework using
#Fermi-LAT Pass 8 gamma-ray data from dwarf spheroidal galaxies.

#It evaluates whether the observed signal amplitude in the
#18–60 GeV range correlates with halo mass (J-factor).

# ------------------------------------------------------------

galaxy_masses = {
    "Leo_2": 17.5,
    "Canes_Venatici_1": 17.2,
    "Ursa_Major_1": 17.2,
    "Ursa_Major_2": 17.9,
    "Canes_Venatici_2": 17.1,
    "Leo_4": 16.3,
    "Hercules": 16.9,
    "Segue_2": 16.2,
    "Triangulum_2": 19.1,
    "Hydra_2": 16.2,
    "Eridanus_2": 15.6,
    "Aquarius_3": 16.0,
    "Segue_1": 19.5,
    "Willman_1": 19.1,
    "Draco": 18.8,
    "Coma_Berenices": 19.0,
    "Reticulum_2": 18.9,
    "Bootes_1": 18.8,
    "Sculptor": 18.6,
    "Fornax": 18.2,
    "Carina": 18.1,
    "Leo_1": 17.7
}

# Energy window where the RIC signal is expected
E_MIN = 18.0  # GeV
E_MAX = 60.0  # GeV


def analyze_mass_invariance():
    """
    Tests whether the gamma-ray signal in the RIC energy window
    is correlated with halo mass (J-factor).

    RIC prediction:
        Signal amplitude should be invariant with respect to mass
        once saturation is reached.

    Standard particle DM expectation:
        Signal should scale with mass / J-factor.
    """

    # Search for FITS files locally or in ./results
    files = glob.glob("*_final.fits")
    if not files:
        files = glob.glob("results/*_final.fits")

    if not files:
        print("No *_final.fits files found.")
        return

    masses = []
    signal_counts = []
    names = []

    print("\n--- RIC MASS INVARIANCE CHECK (Batch 2) ---")
    print(f"Energy window: {E_MIN} – {E_MAX} GeV")

    for f in files:
        basename = os.path.basename(f)
        name_key = basename.replace("_final.fits", "")

        if name_key not in galaxy_masses:
            print(f"Warning: No mass (J-factor) defined for '{name_key}'. Skipping.")
            continue

        mass = galaxy_masses[name_key]

        try:
            with fits.open(f) as hdul:
                # Load energies (MeV → GeV)
                energies = hdul[1].data["ENERGY"] / 1000.0

                # Count photons in the RIC energy window
                mask = (energies >= E_MIN) & (energies <= E_MAX)
                count = np.sum(mask)

                masses.append(mass)
                signal_counts.append(count)
                names.append(name_key)

                print(f"  > {name_key:<20} | J-factor: {mass:.2f} | Signal: {count} photons")

        except Exception as e:
            print(f"Error processing {basename}: {e}")

    if len(masses) < 3:
        print("Not enough data points for statistical analysis.")
        return

    # ------------------------------------------------------------
    # Plot
    # ------------------------------------------------------------

    plt.figure(figsize=(10, 6))

    plt.scatter(
        masses,
        signal_counts,
        c="red",
        s=120,
        edgecolors="black",
        label="dSph observations",
        zorder=3
    )

    # Linear regression through the data (trend visualization only)
    z = np.polyfit(masses, signal_counts, 1)
    p = np.poly1d(z)
    plt.plot(
        masses,
        p(masses),
        "b-",
        linewidth=2,
        alpha=0.7,
        label=f"Linear trend (slope = {z[0]:.2f})"
    )

    for i, name in enumerate(names):
        plt.annotate(
            name,
            (masses[i], signal_counts[i]),
            xytext=(5, 5),
            textcoords="offset points",
            fontsize=8
        )

    plt.title(
        f"RIC Test: Mass Invariance ({len(names)} Galaxies)\n"
        "Is the signal independent of the J-factor?",
        fontsize=12
    )
    plt.xlabel("Galaxy Mass Proxy (log10 J-factor)", fontsize=11)
    plt.ylabel(f"Photon Counts in {E_MIN}–{E_MAX} GeV", fontsize=11)
    plt.grid(True, alpha=0.3)
    plt.legend()

    output_file = "RIC_MASS_INVARIANCE_BATCH2.png"
    plt.savefig(output_file, dpi=300)
    print(f"\nPlot saved as: {output_file}")

    # ------------------------------------------------------------
    # Statistics
    # ------------------------------------------------------------

    corr, p_value = pearsonr(masses, signal_counts)

    print("\nSTATISTICAL RESULTS:")
    print(f"   Pearson correlation coefficient (R): {corr:.3f}")
    print(f"   Linear slope: {z[0]:.3f}")

    print("\nINTERPRETATION:")
    if corr > 0.5:
        print("Strong mass dependence detected (consistent with particle DM scaling).")
    elif abs(corr) < 0.3:
        print("No significant correlation detected.")
        print("   Consistent with RIC saturation / mass invariance.")
    else:
        print("Weak trend detected, not statistically decisive.")


if __name__ == "__main__":
    analyze_mass_invariance()
