"""Example usage with a local file."""

from modis_nd_bias import process_files

files = [
    "/home/iris/Documents/Github/eScience-Tools-in-Climate-Science/"
    "project-01-modis-nd-bias/data/open/MYD06_L2.A2010049.2035.061.2018058173949.hdf"
]
process_files(files, "./processed", resolutions=[1, 5, 20])
