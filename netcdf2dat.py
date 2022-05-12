from argparse import ArgumentParser
from netCDF4 import Dataset
import numpy as np

argparser = ArgumentParser(description="Convert mie netCDF to .dat file.")
argparser.add_argument("input", help="Input file path.")
argparser.add_argument("output", help="Output file path.")
args = argparser.parse_args()

nc_data = Dataset(args.input, "r", format="NETCDF4")

print(nc_data)
print(np.array(nc_data["phase"][:]).shape)

macke = open(args.output, "w")

ntheta = nc_data["ntheta"][0, 0, :].item()
nmom = nc_data["nmom"][0, 0, :].item()
nlam = len(nc_data.dimensions["nlam"])
nreff = len(nc_data.dimensions["nreff"])

macke.write("\n  ")
for r in range(nreff):
    macke.write("%.5g  " % nc_data["reff"][r])
macke.write(" Micron mean \"spheres\" size (A. Macke, 2001)\n")
macke.write("     " + str(ntheta) + "     " + str(nmom - 1) + "      " + str(nlam))
macke.write("\n Angles [DEGREES]\n")
theta = np.array(nc_data["theta"])[0, 0, :, :ntheta].flatten()[::-1].tolist()
i = 0
for t in theta:
    macke.write("%.5g  " % t)
    i += 1
    if (i >= 11):
        i = 0
        macke.write("\n")
macke.write("\n")
for l in range(nlam):
    for r in range(nreff):
        ntheta = nc_data["ntheta"][l, r, :].item()
        nmom = nc_data["nmom"][l, r, :].item()
        lam = np.array(nc_data["wavelen"]).flatten().tolist()[l]
        reff = np.array(nc_data["reff"]).flatten().tolist()[r]
        ext = np.array(nc_data["ext"][l, r]).item()
        ssa = np.array(nc_data["ssa"][l, r]).item()
        macke.write("     " + str(lam) + "    " + str(ext) + "     " + str(ssa) + " !")
        macke.write("\nPhase Function:\n")
        phase = np.array(nc_data["phase"])[l, r, :, :ntheta].flatten()[::-1].tolist()
        i = 0
        for p in phase:
            macke.write("%.5E  " % p)
            i += 1
            if (i >= 6):
                i = 0
                macke.write("\n")
        macke.write("\nLegendre Polynomial Expansion (divided by 2n+1)\n")
        pmom = np.array(nc_data["pmom"])[l, r, :, :nmom][::-1].flatten()
        ran = 2*np.array(range(nmom))+1
        pmom = (pmom/ran).tolist()
        i = 0
        for p in pmom:
            macke.write("%.5g  " % p)
            i += 1
            if (i >= 10):
                i = 0
                macke.write("\n")
        macke.write("\n")
macke.close()
