from argparse import ArgumentParser
from copy import deepcopy
import numpy as np
import numbers
import json
import re

from netCDF4 import Dataset

argparser = ArgumentParser(description="Parse between JSON and netCDF formats.")
argparser.add_argument("input", help="Input file path.")
argparser.add_argument("output", help="Output file path.")
args = argparser.parse_args()

def parse_json_to_netcdf(json, nc_root, hierarchy):
    hierarchy = deepcopy(hierarchy)
    cur_group = nc_root["/" + "/".join(hierarchy)] if len(hierarchy) > 0 else nc_root

    for name, data in json.items():
        if isinstance(data, dict):
            nc_root.createGroup("/" + "/".join(hierarchy + [name]))
            parse_json_to_netcdf(data, nc_root, hierarchy + [name])
        elif isinstance(data, list) and not all(isinstance(x, numbers.Number) for x in data):
            for obj, i in zip(data, range(len(data))):
                nc_root.createGroup("/" + "/".join(hierarchy + [name + "[" + str(i) + "]"]))
                parse_json_to_netcdf(obj, nc_root, hierarchy + [name + "[" + str(i) + "]"])
        elif isinstance(data, list):
            np_data = np.array(data)
            dims = np_data.shape
            dim_names = ["DIM" + str(d) + "_" + name for d in range(len(dims))]
            for dim_name, dim in zip(dim_names, dims):
                cur_group.createDimension(dim_name, dim);
            nc_var = nc_root.createVariable("/" + "/".join(hierarchy + [name]), np_data.dtype, tuple(dim_names))
            nc_var[:] = np_data
        elif isinstance(data, numbers.Number) and not isinstance(data, bool):
            np_data = np.array(data)
            nc_var = nc_root.createVariable("/" + "/".join(hierarchy + [name]), np_data.dtype)
            nc_var[:] = np_data
        else:
            try:
                setattr(cur_group, name, data)
            except:
                setattr(cur_group, name, str(data))

def parse_netcdf_to_json(nc_root, hierarchy):
    hierarchy = deepcopy(hierarchy)
    cur_group = nc_root["/" + "/".join(hierarchy)] if len(hierarchy) > 0 else nc_root
    json = {}
    for attr_name in cur_group.ncattrs():
        json[attr_name] = getattr(cur_group, attr_name)
    for variable in cur_group.variables:
        np_arr = np.array(cur_group[variable][:])
        if np_arr.shape == ():
            json[variable] = np_arr[()].item()
        else:
            json[variable] = np_arr.tolist()
    for group in cur_group.groups:
        sub_json = parse_netcdf_to_json(nc_root, hierarchy + [group])
        try:
            index = int(group.split("[")[1].split("]")[0])
        except:
            index = None
        if index == 0:
            json[group.split("[")[0]] = [sub_json]
        elif not index is None:
            json[group.split("[")[0]] = json[group.split("[")[0]] + [sub_json]
        else:
            json[group] = sub_json
    return json

def walktree(root):
    yield root.groups.values()
    for value in root.groups.values():
        yield from walktree(value)

if args.input.endswith(".json") and args.output.endswith(".nc"):
    with open(args.input) as json_file:
        json_data = json.loads(json_file.read())
    nc_data = Dataset(args.output, "w", format="NETCDF4")
    parse_json_to_netcdf(json_data, nc_data, [])
    for children in walktree(nc_data):
        for child in children:
            print(child)
    nc_data.close()
elif args.input.endswith(".nc") and args.output.endswith(".json"):
    nc_data = Dataset(args.input, "r", format="NETCDF4")
    json_data = parse_netcdf_to_json(nc_data, [])
    nc_data.close()
    with open(args.output, "w") as json_file:
        json.dump(json_data, json_file, indent=4)
else:
    print("Please either specify an input JSON and output netCDF or an input netCDF and output JSON.")
