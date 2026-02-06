# MODIS Nd Bias Package

This package contains the core implementation for the MODIS Nd bias analysis,
plus a CLI that processes MODIS HDF files and generates the same plots as the
original notebook.

## Install (editable)

```bash
cd /home/iris/Documents/Github/eScience-Tools-in-Climate-Science/project-01-modis-nd-bias
python -m pip install -e .[cli,plot]
```

## CLI Usage

Process all HDFs under the default data folders (`data/open`, `data/closed`,
`data/disorganized`) and write NetCDF + plots to `data/processed`:

```bash
python -m modis_nd_bias.cli --recursive --plots
```

Process a specific folder:

```bash
python -m modis_nd_bias.cli --input-dir /path/to/hdf --recursive --plots
```

Process multiple folders:

```bash
python -m modis_nd_bias.cli --input-dirs /path/open /path/closed --recursive --plots
```

Write outputs somewhere else (plots are saved in the same output directory):

```bash
python -m modis_nd_bias.cli --recursive --plots --output-dir /path/to/output
```

Generate plots only (expects `*_processed.nc` files in the output directory):

```bash
python -m modis_nd_bias.cli --plots-only --output-dir /path/to/output
```

## Notes

- Plotting requires `matplotlib` (`pip install -e .[plot]`).
- HDF4 reading requires GDAL (`pip install -e .[cli]`).
- If you don’t install the package, you can run with:

```bash
PYTHONPATH=/home/iris/Documents/Github/eScience-Tools-in-Climate-Science/project-01-modis-nd-bias/src \
python -m modis_nd_bias.cli --recursive --plots
```
