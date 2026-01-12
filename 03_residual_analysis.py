import astropy.io.fits as fits
import matplotlib.pyplot as plt
import numpy as np
import glob
import os
from scipy.optimize import curve_fit

# ------------------------------------------------------------
# Background model: simple power law
# Used to approximate the diffuse gamma-ray background
# ------------------------------------------------------------

def power_law(x, a, k):
    return a * x**(-k)


def analyze_residuals():
    """
    Performs a sideband background fit and computes residuals
    for stacked dSph gamma-ray data.

    The signal window (18–60 GeV) is excluded from the fit.
    Residuals are defined as:
        (Data - Model) / Model

    This test is used to identify structured excesses
    inconsistent with a smooth background.
    """

    # Search for final FITS files in the current directory
    files = glob.glob("*_final.fits")

    # Fallback: allow execution from one directory above
    if not files:
        files = glob.glob("results/*_final.fits")

    if not files:
        print("No *_final.fits files found. Are you in the correct directory?")
        return

    print(f"Computing residuals (sideband fit) for {len(files)} galaxies...")

    # ------------------------------------------------------------
    # Collect photon energies from all targets
    # ------------------------------------------------------------

    all_energies = []

    for f in files:
        try:
            with fits.open(f) as hdul:
                # ENERGY is assumed to be in MeV → convert to GeV
                all_energies.extend(hdul[1].data["ENERGY"] / 1000.0)
        except Exception as e:
            print(f"Error reading {f}: {e}")

    if not all_energies:
        print("No photon data loaded.")
        return

    data = np.array(all_energies)
    n_photons = len(data)
    print(f"Total events: {n_photons}")

    # ------------------------------------------------------------
    # Energy binning (logarithmic)
    # Slightly wider bins to suppress statistical noise
    # ------------------------------------------------------------

    bins = np.logspace(np.log10(10), np.log10(300), 12)
    counts, edges = np.histogram(data, bins=bins)

    bin_centers = np.sqrt(edges[:-1] * edges[1:])
    bin_widths = edges[1:] - edges[:-1]

    # Flux estimate: counts per GeV
    with np.errstate(divide="ignore", invalid="ignore"):
        flux = counts / bin_widths
        errors = np.sqrt(counts) / bin_widths

    # Clean NaN / Inf values (e.g. empty bins)
    flux = np.nan_to_num(flux)
    errors = np.nan_to_num(errors)

    # ------------------------------------------------------------
    # Sideband fit
    # Exclude the signal window (18–60 GeV)
    # ------------------------------------------------------------

    mask_background = (bin_centers < 18.0) | (bin_centers > 60.0)

    if np.sum(mask_background) < 3:
        print("Not enough sideband data points for a stable fit.")
        return

    try:
        popt, pcov = curve_fit(
            power_law,
            bin_centers[mask_background],
            flux[mask_background],
            p0=[1000, 2.5],
            maxfev=10000
        )
        print(f"Background fit successful: spectral index = {popt[1]:.2f}")
    except Exception as e:
        print(f"Background fit failed: {e}")
        return

    # Evaluate background model over full energy range
    model_flux = power_law(bin_centers, *popt)

    # ------------------------------------------------------------
    # Residuals: (Data - Model) / Model
    # ------------------------------------------------------------

    with np.errstate(divide="ignore", invalid="ignore"):
        residuals = (flux - model_flux) / model_flux
        residual_errors = errors / model_flux

    # ------------------------------------------------------------
    # Plot
    # ------------------------------------------------------------

    plt.figure(figsize=(10, 6))

    # Zero line (pure background expectation)
    plt.axhline(0.0, color="black", linewidth=1, alpha=0.5)

    # RIC prediction marker
    plt.axvline(
        20.7,
        color="red",
        linestyle="--",
        alpha=0.8,
        label="RIC prediction (20.7 GeV)"
    )

    # Residual data points
    plt.errorbar(
        bin_centers,
        residuals,
        yerr=residual_errors,
        fmt="o",
        color="blue",
        ecolor="lightblue",
        elinewidth=2,
        capsize=4,
        markersize=6,
        label="Residuals (Data − Model) / Model"
    )

    # Visual indication of the masked signal window
    plt.axvspan(
        18.0,
        60.0,
        color="red",
        alpha=0.05,
        label="Signal window (excluded from fit)"
    )

    plt.title(
        f"RIC Residual Test (N={len(files)} dSphs, {n_photons} photons)\n"
        "Sideband background fit (<18 & >60 GeV)",
        fontsize=12
    )
    plt.xlabel("Energy (GeV)", fontsize=11)
    plt.ylabel("Residuals ((Data − Model) / Model)", fontsize=11)
    plt.xscale("log")
    plt.grid(True, which="both", alpha=0.3)
    plt.legend()

    output_file = "RIC_RESIDUALS_CHECK.png"
    plt.savefig(output_file, dpi=300)
    print(f"\nFigure saved as: {output_file}")


if __name__ == "__main__":
    analyze_residuals()
