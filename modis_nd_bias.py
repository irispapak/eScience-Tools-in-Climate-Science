#!/usr/bin/env python
# coding: utf-8

"""
Investigating biases in droplet number concentration due to averaging of raw satellite observations.
Refactored from Jupyter Notebook for CLI usage.
"""

import glob
import logging
import os
import sys
import warnings
import argparse
from typing import List, Tuple, Optional, Any

import numpy as np
import xarray as xr

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[logging.StreamHandler(sys.stdout)]
)
logger = logging.getLogger(__name__)

# Optional dependencies check
# Optional dependencies check
try:
    from osgeo import gdal
except ImportError:
    logger.warning("Warning: 'osgeo' (GDAL) not found. HDF4 handling might fail.")

try:
    from netCDF4 import Dataset
except ImportError:
    logger.warning("Warning: 'netCDF4' not found. Saving to NetCDF might fail.")

try:
    import statsmodels.api as sm
except ImportError:
    logger.warning("Warning: 'statsmodels' not found. Analysis might be limited.")

# Matplotlib import removed as it is currently unused.

# xr.set_options(display_style='html') # Not needed for CLI
# warnings.filterwarnings('ignore') # Better to handle specific warnings

# Constants
VARNAMES1 = ['Cloud_Fraction', 'Cloud_Top_Temperature']
VARNAMES2 = ['Cloud_Effective_Radius', 'Cloud_Optical_Thickness', 'Cloud_Water_Path', 'Cloud_Multi_Layer_Flag']
VARNAMES3 = ['Cloud_Mask_1km']
DEFAULT_RESOLUTIONS = [1, 5, 10, 20, 25, 50, 100]


MISSING_FLOAT = -9999.0
MISSING_INT = -32768
BAD_VALUE = -999.0


def saturation_vapor_pressure(T: np.ndarray) -> np.ndarray:
    """
    Calculates the saturation vapor pressure (es) in Pascals.
    
    Examples:
        >>> np.isclose(saturation_vapor_pressure(273.15), 611.2, atol=1.0)
        True
        >>> np.isclose(saturation_vapor_pressure(373.15), 101325, rtol=0.05) # Approx boiling point
        True
    """
    # formula from Bolton (1980) or similar?
    # Keeping original formula structure
    return np.exp(54.84 - 6763.22/T - 4.21*np.log(T) + 0.00037*T + np.tanh(0.0415*(T-218.8)) *(53.878 - 1331.22/T - 9.44523*np.log(T) + 0.014025*T))


def saturation_adiabatic_lapse_rate(T: np.ndarray, p: np.ndarray) -> np.ndarray:
    """
    Calculates the moist adiabatic lapse rate (Gamma_s).
    
    Examples:
        >>> # Lapse rate should be positive and less than dry adiabatic (~9.8 K/km)
        >>> gamma = saturation_adiabatic_lapse_rate(288.15, 101325)
        >>> 0.0 < gamma < 0.0098
        True
    """
    g = 9.81 
    H = 2501000  # Latent heat of vaporization approx?
    ep = 0.622   # Ratio of gas constants
    Rsd = 287    # Specific gas constant for dry air
    cpd = 1003.5 # Specific heat of dry air
    
    es_val = saturation_vapor_pressure(T)
    numerator = g * (1 + ((H * ep * es_val) / (Rsd * (p - es_val) * T)))
    denominator = cpd + ((H ** 2 * ep ** 2 * es_val) / (Rsd * (p - es_val) * T ** 2))
    return numerator / denominator


def calculate_droplet_concentration(re: xr.DataArray, tau: xr.DataArray, T: float = 275.0, P: float = 95000.0, f_ad: float = 0.8) -> xr.DataArray:
    """
    Calculates the droplet number concentration (Nd), given the cloud optical
    thickness (tau) and cloud effective radius (re).
    
    Args:
        re: Effective radius (microns)
        tau: Optical thickness
        T: Temperature (K)
        P: Pressure (Pa)
        f_ad: Adiabatic fraction
        
    Returns:
        Nd: Droplet number concentration
    """
    Qext = 2            
    ro_w = 997 * 10**3     # Water density [g/m^3]

    g = 9.81
    Cp = 1004              # [J/kg K]
    ro_a = 1.2             # Air density [kg/m^3]
    Lv =  2.5 * 10**6      # Latent heat of vaporization [J/kg]
    gamma_d = g/Cp
    
    lapse_rate = saturation_adiabatic_lapse_rate(np.array(T), np.array(P))
    Cw = f_ad * ((ro_a * Cp * (gamma_d - lapse_rate) / (Lv)) * 1000)  # [g/m^4]
    
    gamma = ((5**0.5)/(2*np.pi*0.8)) * (Cw/(Qext*ro_w))**0.5

    # Safe handling for re being 0 or nan
    N = (gamma * tau ** 0.5 * (re * (1e-6)) ** (-5. / 2)) * 1e-6
    return N


def get_subdataset(file: str, varnames: List[str]) -> List[str]:
    """Retrieves subdataset names from an HDF file using GDAL."""
    g = gdal.Open(file)
    if g is None:
        raise IOError(f"Could not open file with GDAL: {file}")
    subdatasets = g.GetSubDatasets()
    l = []
    for varname in varnames:
        # Check against the identifier string in the tuple (name, description)
        l.extend([s[0] for s in subdatasets if varname in s[0].split(':')[-1]])
    return l


def read_hdf4(file: str, varnames: List[str]) -> xr.Dataset:
    """Reads specific variables from an HDF4 file into an xarray Dataset."""
    ld = []
    fname_list = get_subdataset(file, varnames)

    for fname in fname_list:
        # Note: open_rasterio is deprecated, but currently required for this GDAL-based workflow.
        # Consider migrating to rioxarray in the future.
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", category=DeprecationWarning)
            # Use chunks=None or explicit chunks to possibly return dask array if needed, 
            # but current logic assumes in-memory.
            myDataset = xr.open_rasterio(fname)
            
        short_name = fname.split(':')[-1]
        ld.append(myDataset.to_dataset(name=short_name))

    return xr.merge(ld)


def coarsen_image(da, avg_window_size=10, boundary='exact'):
    """Coarsens an xarray.DataArray based on the averaging window size."""
    return da.coarsen(x=avg_window_size, y=avg_window_size, boundary=boundary).mean(skipna=True)


def slice_image(da, x0, x1, y0, y1):
    """Slices an xarray.DataArray based on Along_Swath and Across_Swath"""
    return da.where((((da.y >= y0) & (da.y <= y1)) & 
                        ((da.x >= x0) & (da.x <= x1))), drop=True)


def store_cloud_data(filename, var1, var2, var3):
    """Reads the hdf file and extracts cloud information."""
    ds1 = read_hdf4(filename, var1)
    ds2 = read_hdf4(filename, var2)
    ds3 = read_hdf4(filename, var3)

    dim_x_fine = ds2['x'].values
    dim_y_fine = ds2['y'].values
    
    # store variables
    cf = ds1['Cloud_Fraction']
    top_temp = ds1['Cloud_Top_Temperature']
    reff = ds2['Cloud_Effective_Radius']
    tau = ds2['Cloud_Optical_Thickness']
    lwp = ds2['Cloud_Water_Path']
    multi_layer = ds2['Cloud_Multi_Layer_Flag']
    cld_mask_init = ds3['Cloud_Mask_1km'][:,:,0]

    # choose only marine clouds 
    cld_mask_copy = np.where((cld_mask_init==57) | (cld_mask_init==41), 1, 0)

    # add 1 extra dimension to match the other variables
    cld_mask_copy_newdim = np.expand_dims(cld_mask_copy, axis=0)

    # create dataarray with correct dimensions
    cld_mask = xr.DataArray(name='Cloud_Mask_1km', data=cld_mask_copy_newdim, dims=['band','y','x'], coords=[ [1], dim_y_fine, dim_x_fine ] )

    del cld_mask_copy, cld_mask_copy_newdim, cld_mask_init

    return dim_x_fine, dim_y_fine, cf, top_temp, reff, tau, lwp, multi_layer, cld_mask


def ignore_missing(top_temp: xr.DataArray, reff: xr.DataArray, tau: xr.DataArray, lwp: xr.DataArray) -> Tuple[xr.DataArray, xr.DataArray, xr.DataArray, xr.DataArray]:
    """Set all missing values (e.g., -9999) to NaNs."""
    # Use constants defined at module level
    top_temp = top_temp.where(top_temp != MISSING_INT)
    reff = reff.where(reff != MISSING_FLOAT)
    tau = tau.where(tau != MISSING_FLOAT)
    lwp = lwp.where(lwp != MISSING_FLOAT)

    return top_temp, reff, tau, lwp


def interpolation_from_5km_to_1km(dimx, dimy, cf, t):
    """Changes coordinates of 5km grid to 1km grid and interpolates."""
    cf_coords = cf.assign_coords(y = (cf.y-0.5)*5, x = (cf.x-0.5)*5) 
    t_coords = t.assign_coords(y = (t.y-0.5)*5, x = (t.x-0.5)*5) 

    cf_interp = cf_coords.interp(x = dimx, y = dimy) 
    t_interp = t_coords.interp(x = dimx, y = dimy) 

    del cf, t, cf_coords, t_coords
    return cf_interp, t_interp


def scaling_variables(reff, tau, cf_interp, top_temp_interp):
    """Scales the variables according to HDF4 metadata."""
    scale_factor = 0.00999999977648258
    add_offset = -15000

    reff = reff*scale_factor
    tau = tau*scale_factor
    cf = cf_interp*scale_factor
    top_temp = (top_temp_interp-add_offset)*scale_factor

    del top_temp_interp, cf_interp
    return reff, tau, cf, top_temp


def filtering(cld_mask, tau, reff, lwp, cf, top_temp, multi_layer):
    """Filters values for warm, single-layer marine clouds with sufficient optical thickness."""
    # tau>5, cld_mask==1 (marine), top_temp>273 (warm), reff>5um, multi_layer=1 (single)
    gen_mask = xr.where( (tau>5) & (cld_mask==1) & (top_temp>273) & (reff>5) & (multi_layer==1), cld_mask, 0)

    tau = xr.where(gen_mask==1, tau, np.nan)
    reff = xr.where(gen_mask==1, reff, np.nan)
    lwp = xr.where(gen_mask==1, lwp, np.nan)
    # top_temp isn't returned but used for mask generation
    
    del top_temp, multi_layer, cld_mask, gen_mask
    return tau, reff, lwp, cf


def divide_in_100X100km_regions(dimx, dimy, reff, tau, lwp, cf):
    """Cuts image into 100x100px regions with high cloud fraction."""
    reg = 100
    if len(dimx) < reg or len(dimy) < reg:
        return [], [], [], [], 0

    int_dim_x = int(len(dimx)/reg)*reg
    int_dim_y = int(len(dimy)/reg)*reg

    reff_reg = []
    tau_reg = []
    nd_reg = []
    lwp_reg = []

    ct = 0
    # Step size is usually just 'reg', the original code had 'reg+1' which assumes specific skipping? 
    # Keeping original logic for fidelity unless it causes issues.
    step = reg + 1 
    
    for i in range(0, int_dim_x, step):
        for j in range(0, int_dim_y, step):
            # Ensure slicing bounds don't exceed dimensions
            # Original code used i+2, i+2+reg. Assuming margins.
            
            # Simple bounds check
            if i+2+reg > len(dimx) or j+2+reg > len(dimy):
                continue

            cf_now = slice_image(cf, i+2, i+2+reg, j+2, j+2+reg) 
            cf_regmean = cf_now.mean(skipna=True)

            if cf_regmean >= 0.99:
                reff_now = slice_image(reff, i+2, i+2+reg, j+2, j+2+reg)
                tau_now = slice_image(tau, i+2, i+2+reg, j+2, j+2+reg)
                lwp_now = slice_image(lwp, i+2, i+2+reg, j+2, j+2+reg)

                nd_now = calculate_droplet_concentration(reff_now, tau_now)

                reff_reg.append(reff_now)
                tau_reg.append(tau_now)
                lwp_reg.append(lwp_now)
                nd_reg.append(nd_now)
                ct += 1

    del tau, reff, lwp, cf
    return reff_reg, tau_reg, nd_reg, lwp_reg, ct


def coarsening_and_regional_mean(resol, nresol, nreg, nd_reg_fine, lwp_reg_fine, reff_reg_fine, tau_reg_fine):
    """Coarsens resolution and calculates regional means."""
    nd_coarse_regmean = np.nan*np.zeros((nreg, nresol))
    lwp_coarse_regmean = np.nan*np.zeros((nreg, nresol))
    lwp_coarse_regstd = np.nan*np.zeros((nreg, nresol))
    reff_coarse_regmean = np.nan*np.zeros((nreg, nresol))

    for r in range(nreg):
        for i in range(nresol):
            if i == 0: # 1km resolution (finest)
                lwp_coarse_regmean[r,0] = lwp_reg_fine[r].mean(skipna=True)
                lwp_coarse_regstd[r,0] = lwp_reg_fine[r].std(skipna=True)
                nd_coarse_regmean[r,0] = nd_reg_fine[r].mean(skipna=True)
                reff_coarse_regmean[r,0] = reff_reg_fine[r].mean(skipna=True)
            else:
                # Coarsen
                lwp_reg_coarse = coarsen_image(lwp_reg_fine[r], resol[i])
                reff_reg_coarse = coarsen_image(reff_reg_fine[r], resol[i])
                tau_reg_coarse = coarsen_image(tau_reg_fine[r], resol[i])
                nd_reg_coarse = calculate_droplet_concentration(reff_reg_coarse, tau_reg_coarse)

                lwp_coarse_regmean[r,i] = lwp_reg_coarse.mean(skipna=True)
                lwp_coarse_regstd[r,i] = lwp_reg_coarse.std(skipna=True)
                nd_coarse_regmean[r,i] = nd_reg_coarse.mean(skipna=True)
                reff_coarse_regmean[r,i] = reff_reg_coarse.mean(skipna=True)

    return nd_coarse_regmean, lwp_coarse_regmean, lwp_coarse_regstd, reff_coarse_regmean


def calculate_bias_nd_and_homogeneity(nresol, nreg, nd_regmean, lwp_regmean, lwp_regstd):
    """Calculates bias in Nd and homogeneity measure."""
    bias = np.nan*np.zeros((nreg, nresol))
    homog = np.nan*np.zeros((nreg))

    for r in range(nreg):
        for i in range(nresol):
            if nd_regmean[r,0] != 0 and not np.isnan(nd_regmean[r,0]):
                bias[r,i] = (nd_regmean[r,i]-nd_regmean[r,0])*100/nd_regmean[r,0]
        
        if lwp_regmean[r,0] != 0 and not np.isnan(lwp_regmean[r,0]):
            homog[r] = lwp_regstd[r,0]/lwp_regmean[r,0]    

    return bias, homog


def save_info_to_nc(output_dir, file, resol, nresol, nreg, bias, homog, nd_regmean, lwp_regmean, reff_regmean):
    """Saves processed data to NetCDF."""
    input_filename = os.path.basename(file)
    filename_base = os.path.splitext(input_filename)[0]
    output_path = os.path.join(output_dir, filename_base + '_processed.nc')

    logger.info(f"Saving to {output_path}...")
    ncout = Dataset(output_path, 'w', format='NETCDF4')

    # create dimensions
    ncout.createDimension('reg', nreg)
    ncout.createDimension('resol', nresol)

    # create variables
    region = ncout.createVariable('reg', np.int64, ('reg',))
    resolution = ncout.createVariable('resol', np.int64, ('resol',))

    bias_nd = ncout.createVariable('bias_nd', np.float32, ('reg','resol'))
    homogeneity = ncout.createVariable('homog', np.float32, ('reg',))
    nd_mean = ncout.createVariable('nd_mean', np.float32, ('reg','resol'))
    lwp_mean = ncout.createVariable('lwp_mean', np.float32, ('reg', 'resol'))
    reff_mean = ncout.createVariable('reff_mean', np.float32, ('reg', 'resol'))

    # Attributes
    region.units = '-'
    resolution.units = 'km'
    bias_nd.units = '%'
    homogeneity.units = '-'
    nd_mean.units = '#/cm^3'
    lwp_mean.units = 'g/m^2'
    reff_mean.units = 'um'

    region[:] = np.arange(1, nreg+1, 1)
    resolution[:] = resol[:]
    bias_nd[:,:] = bias[:,:]
    homogeneity[:] = homog[:]
    nd_mean[:,:] = nd_regmean[:,:]
    lwp_mean[:,:] = lwp_regmean[:,:]
    reff_mean[:,:] = reff_regmean[:,:]

    ncout.close()


def process_files(input_files: List[str], output_dir: str, resolutions: List[int]):
    """Main processing loop."""
    if not input_files:
        logger.warning("No HDF files found.")
        return

    nresolutions = len(resolutions)
    total_num_regions = 0

    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    for f in input_files:
        logger.info(f"Processing {f}...")
        try:
            dimx, dimy, cf, top_temp, reff, tau, lwp, multi_layer, cld_mask = store_cloud_data(f, VARNAMES1, VARNAMES2, VARNAMES3)
            
            top_temp, reff, tau, lwp = ignore_missing(top_temp, reff, tau, lwp)
            cf_int, top_temp_int = interpolation_from_5km_to_1km(dimx, dimy, cf, top_temp)
            reff, tau, cf, top_temp = scaling_variables(reff, tau, cf_int, top_temp_int)
            tau, reff, lwp, cf = filtering(cld_mask, tau, reff, lwp, cf, top_temp, multi_layer)

            reff_1km_reg, tau_1km_reg, nd_1km_reg, lwp_1km_reg, nregions = divide_in_100X100km_regions(dimx, dimy, reff, tau, lwp, cf)

            if nregions == 0:
                logger.info(f"  No valid regions found in {os.path.basename(f)}")
                continue

            logger.info(f"  Found {nregions} valid regions.")
            total_num_regions += nregions

            nd_all_regmean, lwp_all_regmean, lwp_all_regstd, reff_all_regmean = coarsening_and_regional_mean(
                resolutions, nresolutions, nregions, nd_1km_reg, lwp_1km_reg, reff_1km_reg, tau_1km_reg
            )

            bias_nd, meas_homog = calculate_bias_nd_and_homogeneity(
                nresolutions, nregions, nd_all_regmean, lwp_all_regmean, lwp_all_regstd
            )

            save_info_to_nc(
                output_dir, f, resolutions, nresolutions, nregions, bias_nd, meas_homog, 
                nd_all_regmean, lwp_all_regmean, reff_all_regmean
            )

        except Exception as e:
            logger.error(f"Error processing {f}: {e}", exc_info=True)


def main():
    parser = argparse.ArgumentParser(description="Process MODIS cloud data to analyze Nd bias.")
    parser.add_argument('--input-dir', type=str, required=True, help='Directory containing HDF files (can use recursive globbing).')
    parser.add_argument('--output-dir', type=str, default='./processed', help='Directory to save processed NetCDF files.')
    parser.add_argument('--recursive', action='store_true', help='Search recursively for HDF files in input directory.')
    
    args = parser.parse_args()

    search_pattern = "**/*.hdf" if args.recursive else "*.hdf"
    
    # Check if input directory exists
    if not os.path.exists(args.input_dir):
        logger.error(f"Input directory not found: {args.input_dir}")
        sys.exit(1)

    fileset = [file for file in glob.glob(os.path.join(args.input_dir, search_pattern), recursive=args.recursive)]
    
    logger.info(f"Found {len(fileset)} files in {args.input_dir}")
    
    process_files(fileset, args.output_dir, DEFAULT_RESOLUTIONS)


if __name__ == "__main__":
    main()
