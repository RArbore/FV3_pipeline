import matplotlib.pyplot as plt
from netCDF4 import Dataset
import numpy as np
import re

argparser = ArgumentParser(description="Create graph of phase function output vs angle (operates on output of libRadtran's mie, not FV3 pipeline).")
argparser.add_argument("input", help="Input netCDF file (output from libRadtran's mie).")
argparser.add_argument("output", help="Output directory to write plots into.")
args = argparser.parse_args()

nc_data = Dataset(args.input, "r")

for l in range(len(nc_data.dimensions["nlam"])):
    for r in range(len(nc_data.dimensions["nreff"])):
        phase_nc_data = nc_data["phase"][l, r, 0, :nc_data["ntheta"][l, r, 0]]
        plt.yscale('log')
        plt.plot(nc_data["theta"][l, r, 0, :nc_data["ntheta"][l, r, 0]], phase_nc_data)
        plt.savefig(args.output + "/wavelen_num_" + str(l) + "_eff_rad_num_" + str(r) + ".png")
        plt.clf()
