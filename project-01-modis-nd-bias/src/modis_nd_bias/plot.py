"""Plotting utilities to reproduce notebook figures from processed NetCDF files."""

from __future__ import annotations

import glob
from pathlib import Path
from typing import Tuple

import numpy as np
from netCDF4 import Dataset
import matplotlib.pyplot as plt


def _load_processed(processed_dir: Path) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    files_nc = sorted(glob.glob(str(processed_dir / "*.nc")))
    if not files_nc:
        raise FileNotFoundError(f"No processed NetCDF files found in {processed_dir}")

    nd_means = None
    lwp_means = None
    reff_means = None
    biases = None
    homogeneities = None

    for i, f in enumerate(files_nc):
        data = Dataset(f, "r")
        data_nd = np.array(data["nd_mean"])
        data_lwp = np.array(data["lwp_mean"])
        data_bias = np.array(data["bias_nd"])
        data_homo = np.array(data["homog"])
        data_reff = np.array(data["reff_mean"])
        data.close()

        if i == 0:
            nd_means = data_nd
            lwp_means = data_lwp
            reff_means = data_reff
            biases = data_bias
            homogeneities = data_homo
        else:
            nd_means = np.concatenate((data_nd, nd_means), axis=0)
            lwp_means = np.concatenate((data_lwp, lwp_means), axis=0)
            reff_means = np.concatenate((data_reff, reff_means), axis=0)
            biases = np.concatenate((data_bias, biases), axis=0)
            homogeneities = np.concatenate((data_homo, homogeneities), axis=0)

    return biases, homogeneities, nd_means, lwp_means, reff_means


def _set_plot_style() -> None:
    plt.rc("ytick", labelsize=15)
    plt.rc("xtick", labelsize=15)
    plt.rc("axes", labelsize=20)
    plt.rc("axes", titlesize=20)
    plt.rc("font", size=20)


def plot_bias_boxplots(biases: np.ndarray, nresolutions: int, output_dir: Path) -> Path:
    data = []
    for i in range(1, nresolutions):
        filtered_data = biases[:, i][~np.isnan(biases[:, i])]
        data.append(filtered_data)

    total_num_regions = len(data[0]) if data else 0
    fig, ax = plt.subplots(figsize=(14, 12))
    plt.boxplot(data, showfliers=False, showmeans=True)
    ax.set_xticklabels(["5", "10", "20", "25", "50", "100"])
    plt.ylabel("Relative Bias of $N_d$ (%)")
    plt.xlabel("Resolution (km)")
    plt.title(f"$N_d$-bias with resolution ({total_num_regions} regions)")

    out = output_dir / "bias-nd-boxplots.jpg"
    fig.savefig(out)
    plt.close(fig)
    return out


def plot_bias_vs_homogeneity(
    biases: np.ndarray,
    homogeneities: np.ndarray,
    color_values: np.ndarray,
    resolutions: list[int],
    output_path: Path,
    color_label: str,
    vmin: float | None = None,
    vmax: float | None = None,
) -> Path:
    nresolutions = len(resolutions)
    fig, axes = plt.subplots(nrows=3, ncols=2, sharex=True, sharey=True, figsize=(15, 15))

    row = [0, 0, 1, 1, 2, 2]
    col = [0, 1, 0, 1, 0, 1]

    for i in range(1, nresolutions):
        axes[row[i - 1], col[i - 1]].set_title(f"{resolutions[i]} km", fontsize=20)
        plot = axes[row[i - 1], col[i - 1]].scatter(
            homogeneities[:],
            biases[:, i],
            s=50,
            c=color_values[:, i],
            vmin=vmin,
            vmax=vmax,
        )

        x = homogeneities[:]
        y = biases[:, i]
        mask = np.isfinite(x) & np.isfinite(y)
        if mask.sum() >= 2:
            slope, intercept = np.polyfit(x[mask], y[mask], 1)
            axes[row[i - 1], col[i - 1]].plot(x, slope * x + intercept)

    plt.xticks(np.arange(0, 1.3, 0.2))
    plt.xlim([0, 1.2])
    cb = fig.colorbar(plot, ax=axes[:, :])
    cb.set_label(color_label)
    fig.text(0.43, 0.06, "stdLWP/meanLWP", ha="center")
    fig.text(0.06, 0.5, "Relative Bias of $N_d$ (%)", va="center", rotation="vertical")

    fig.savefig(output_path)
    plt.close(fig)
    return output_path


def generate_all_plots(processed_dir: Path) -> None:
    processed_dir = Path(processed_dir)
    _set_plot_style()

    biases, homogeneities, nd_means, lwp_means, reff_means = _load_processed(processed_dir)
    resolutions = [1, 5, 10, 20, 25, 50, 100]
    nresolutions = len(resolutions)

    plot_bias_boxplots(biases, nresolutions, processed_dir)

    plot_bias_vs_homogeneity(
        biases,
        homogeneities,
        nd_means,
        resolutions,
        processed_dir / "bias-nd-with-homogeneity-cb-nd.jpg",
        "mean Nd (#/cm$^3$)",
        vmin=0,
        vmax=200,
    )

    plot_bias_vs_homogeneity(
        biases,
        homogeneities,
        reff_means,
        resolutions,
        processed_dir / "bias-nd-with-homogeneity-cb-reff.jpg",
        "mean reff (um)",
    )

    plot_bias_vs_homogeneity(
        biases,
        homogeneities,
        lwp_means,
        resolutions,
        processed_dir / "bias-nd-with-homogeneity-cb-lwp.jpg",
        "mean LWP (g/m$^2$)",
        vmin=0,
        vmax=150,
    )
