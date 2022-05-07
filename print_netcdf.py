from argparse import ArgumentParser
from netCDF4 import Dataset

def walktree(top):
    yield top.groups.values()
    for value in top.groups.values():
        yield from walktree(value)

argparser = ArgumentParser(description="Quickly print netCDF file info.")
argparser.add_argument("input", help="Input file path.")
argparser.add_argument("--var", nargs="?", help="Variable to print contents of.")
args = argparser.parse_args()

nc_data = Dataset(args.input, "r", format="NETCDF4")
print("Top level:")
print(nc_data)
if not args.var is None:
    print("Variable " + args.var + ":")
    print(nc_data[args.var][:])
print("")
print("Walk groups:")
for children in walktree(nc_data):
    for child in children:
        print(child)
nc_data.close()
