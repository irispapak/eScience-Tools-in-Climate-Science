# MODIS-Nd-Bias: Satellite Cloud Droplet Analysis

A Python-based tool for analyzing biases in Cloud Droplet Number Concentration ($N_d$) from MODIS satellite data. 

## Features

-   **Robust Data Processing**: Filters and processes MODIS HDF4 Level 2 cloud products.
-   **Physics-Based Derivation**: Calculates $N_d$ using adiabatic assumptions and observation-specific lapse rates.
-   **Aggregated Statistics**: Computes regional means and biases across multiple spatial resolutions (1km to 100km).
-   **Dual Interface**: Use as a standalone Command Line Interface (CLI) or import as a library for Jupyter Notebooks.

## Installation

Ensure you have the following dependencies installed (via `conda` or `pip`):
-   `numpy`
-   `xarray`
-   `gdal` (for HDF4 handling)
-   `netCDF4` (for output generation)

## Usage

### 1. Command Line Interface (CLI)
You can run the script directly from the terminal to process an entire directory of HDF files.

```bash
python3 modis_nd_bias.py --input-dir ./data/ --output-dir ./results/ --recursive
```

**Arguments:**
-   `--input-dir`: Path to directory containing MODIS `.hdf` files.
-   `--output-dir`: (Optional) Directory to save processed `.nc` files. Defaults to `./processed`.
-   `--recursive`: (Optional) Search for files recursively in subdirectories.

### 2. Jupyter Notebook / Library
For interactive analysis, import the module in your Python scripts or Notebooks. See `modis_nd_bias.ipynb` for a complete example.

```python
import modis_nd_bias

# Define paths and run
files = ["file1.hdf", "file2.hdf"]
modis_nd_bias.process_files(files, "./output_dir", resolutions=[1, 5, 20])
```

## Project Structure
-   `modis_nd_bias.py`: Core logic and CLI entry point.
-   `modis_nd_bias.ipynb`: Demonstration notebook.

---
*Based on the eScience Tools in Climate Science winter course (2020).*
