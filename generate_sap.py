from argparse import ArgumentParser
from netCDF4 import Dataset
import numpy as np
import json

argparser = ArgumentParser(description="Generate SAP input for MODTRAN.")
argparser.add_argument("mie_netCDF", help="netCDF output from libRadtran's mie.")
argparser.add_argument("input_json", help="Input JSON file for MODTRAN.")
argparser.add_argument("output", help="Output SAP file path.")
args = argparser.parse_args()

INPUT_NETCDF = args.mie_netCDF
INPUT_JSON = args.input_json
OUTPUT_SAP = args.output

MODTRAN_LAYDIM = 502

json_data = json.loads(open(INPUT_JSON, "r").read());

ZPCLD = json_data["MODTRAN"][0]["MODTRANINPUT"]["AEROSOLS"]["CLDALT"]["ZPCLD"]
CLD = json_data["MODTRAN"][0]["MODTRANINPUT"]["AEROSOLS"]["CLDALT"]["CLD"]
alts = json_data["MODTRAN"][0]["MODTRANINPUT"]["ATMOSPHERE"]["PROFILES"][0]["PROFILE"]

nc_data = Dataset(INPUT_NETCDF, "r", format="NETCDF4");

NWVSAP = len(nc_data.dimensions["nlam"])
NLGSAP = min(nc_data["nmom"][0, 0, 0], MODTRAN_LAYDIM - 1)
NANSAP = nc_data["ntheta"][0, 0, 0]

wavelen = np.array(nc_data["wavelen"][:])
theta = np.array(nc_data["theta"][:])[0, 0, 0, :NANSAP][::-1]
phase = np.array(nc_data["phase"][:])[:, 0, 0, :NANSAP][:, ::-1]
ext = np.array(nc_data["ext"][:])[:, 0]
ssa = np.array(nc_data["ssa"][:])[:, 0]
pmom = np.array(nc_data["pmom"][:])[:, 0, 0, :NLGSAP]
pmom /= 1 + 2 * np.arange(pmom.shape[-1])

f = open(OUTPUT_SAP, "w")

f.write(str(NWVSAP) + " " + str(NLGSAP - 1) + " " + str(NANSAP) + "\n")
f.write("{:.10f}".format(theta[0]))
for t in theta[1:]:
	f.write(" " + "{:.10f}".format(t))
f.write("\n")

cur_cld_i = 0
last_cld = 0.0
num_alts = min(len(alts), MODTRAN_LAYDIM - 1)
for alt in range(num_alts):
	while alts[alt] > ZPCLD[cur_cld_i] and cur_cld_i < len(ZPCLD) - 1:
		cur_cld_i += 1
		last_cld = CLD[cur_cld_i]
	not_zero = alt < num_alts - 1 and last_cld > 0.0
	for lam in range(len(wavelen)):
		bext = ext[lam] if not_zero else 0.0
		bssa = ssa[lam] if not_zero else 0.0
		f.write("{:.10f}".format(alts[alt] * 1000) + " " + "{:.10f}".format(wavelen[lam]) + " " + "{:.10f}".format(bext) + " " + "{:.10f}".format(bssa))
		for legendre in range(len(pmom[lam])):
			f.write(" " + "{:.10f}".format(pmom[lam][legendre] if not_zero else (1.0 if (alt < num_alts - 1 and legendre == 0) else 0.0)))
		for pha in phase[lam]:
			f.write(" " + "{:.10f}".format(pha if alt < num_alts - 1 else 0.0))
		f.write("\n")
