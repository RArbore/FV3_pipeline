from argparse import ArgumentParser
import matplotlib.pyplot as plt
from netCDF4 import Dataset
import numpy as np

argparser = ArgumentParser(description="Create graph of extinction coefficients vs wavelength (operates on output of libRadtran's mie, not FV3 pipeline).")
argparser.add_argument("input", help="Input netCDF file (output from libRadtran's mie).")
argparser.add_argument("output", help="Output directory to write plots into.")
args = argparser.parse_args()

nc_data = Dataset(args.input, "r")
nlam = len(nc_data.dimensions["nlam"])
nreff = len(nc_data.dimensions["nreff"])

for r in range(nreff):
    ext = nc_data["ext"][:, r]
    wavelen = nc_data["wavelen"][:]
    indices = np.argsort(wavelen)

    plt.plot(wavelen[indices], ext[indices])
    plt.xlabel("Wavelength (microns)")
    plt.ylabel("Ext coeff (km-1 m3/g)")
    plt.savefig(args.output + "/eff_rad_num_" + str(r) + ".png")
    plt.clf()
