from argparse import ArgumentParser
from netCDF4 import Dataset
import numpy as np

def diff_groups(root1, root2, path):
	path = "root" if path == "" else path
	
	for k1, k2 in zip(root1.dimensions, root2.dimensions):
		v1 = root1.dimensions[k1]
		v2 = root2.dimensions[k2]
		if not (k1, len(v1)) == (k2, len(v2)):
			print("Dimensions differ: names are " + k1 + " and " + k2 + ", and lengths are " + str(len(v1)) + " and " + str(len(v2)))

	for k1, k2 in zip(root1.variables, root2.variables):
		v1 = root1.variables[k1]
		v2 = root2.variables[k2]
		if not k1 == k2 or not (v1[:] == v2[:]).all():
			diff_sum = np.sum(np.abs(v1[:] - v2[:]))
			print("Variables differ: names are " + k1 + " and " + k2 + ", and the sum of absolute differences is " + str(diff_sum))

	if not (len(root1.groups) == len(root2.groups)):
		print("Different number of groups of " + path)
	for (k1, v1), (k2, v2) in zip(root1.groups, root2.groups):
		if k1 == k2:
			diff_groups(root1[k1], root[k2])
		else:
			print("Differing group names: " + k1 + " and " + k2)

argparser = ArgumentParser(description="Compare 2 netCDF files.")
argparser.add_argument("input1", help="First input file path.")
argparser.add_argument("input2", help="Second input file path.")
args = argparser.parse_args()

nc_data1 = Dataset(args.input1, "r", format="NETCDF4")
nc_data2 = Dataset(args.input2, "r", format="NETCDF4")

diff_groups(nc_data1, nc_data2, "")
