from argparse import ArgumentParser
import matplotlib.pyplot as plt
from netCDF4 import Dataset
import numpy as np
import os

argparser = ArgumentParser(description="Create graph of toa_radiance vs azimuth angle.")
argparser.add_argument("input", help="Input netCDF file (output from FV3 pipeline).")
argparser.add_argument("case", help="Profile to construct plots for.")
argparser.add_argument("sza", help="SZA to construct plots for.")
argparser.add_argument("lam", help="Wavelength to construct plots for.")
argparser.add_argument("output", help="Output plot.")
args = argparser.parse_args()

netcdf = Dataset(args.input, "r", format="NETCDF4")
angles = netcdf.variables["azinp"][args.case, args.sza, :]
lam_pos = netcdf.variables["toa_radiances_wavelengths"][args.case, args.sza, :, :]
lam_pos = np.where(lam_pos == float(args.lam))
toa_rads = netcdf.variables["toa_radiances"][args.case, args.sza, :, :][lam_pos]

plt.clf()
plt.plot(angles, toa_rads)
plt.savefig(args.output)
