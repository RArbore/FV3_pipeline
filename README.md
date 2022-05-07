# FV3 Pipeline
The FV3 pipeline consists of several tools to perform file conversions and MODTRAN calculations.

## Project Files
### FV3_pipeline.sh
Script for running the entire FV3 pipeline on an input netCDF file. This script accepts a single argument: the input netCDF file.
```
./FV3_pipeline.sh <input file>
```
#### Expected format for input netCDF
See the "Input netCDF" section in FORMAT.md for this information

#### Output from FV3 pipeline
Running this script will create intermediary files, as well as the final pipeline output. Each file is prefixed with the input file name, sans the ending ".nc".
The following files are created:
```
- <input file>.json               # conversion between input netCDF and JSON file formats
- <input file>_modtran_input.json # restructure input JSON and add info for MODTRAN
- $SCRATCH/output/*               # the raw output from MODTRAN; stored in your $SCRATCH directory as this output is quite large, and will likely exceed the quota for your home directory
- <input file>_modtran_output.nc  # the final output netCDF; relevant fields are extracted from the raw MODTRAN output and placed in this netCDF
```

#### Final netCDF format
See the "Output netCDF" section in FORMAT.md for this information

### json_netcdf.py
This tool is used for converting between the JSON and netCDF file formats. Although it is only used for netCDF to JSON conversion in the standard pipeline, it also allows for conversion from JSON files to netCDF files.
```
usage: python3 json_netcdf.py [-h] input output

Parse between JSON and netCDF formats.

positional arguments:
  input       Input file path.
  output      Output file path.

optional arguments:
  -h, --help  show this help message and exit
```

### FV3_construct_input.py
This tool is used to construct the input MODTRAN json from the JSON converted directly from the input netCDF. This tool produces the JSON input that is fed directly into MODTRAN.
```
usage: python3 FV3_construct_input.py [-h] input output

Add structure to JSON files for MODTRAN input (specific to FV3 OSSE inputs).

positional arguments:
  input       Input file path.
  output      Output file path.

optional arguments:
  -h, --help  show this help message and exit
```

### mpi_modtran.cpp & Makefile
This tool is what actually calls MODTRAN on input cases. It uses MPI to run multiple cases in parallel.
Due to licensing issues, *you should only run this tool on a single node*. On NERSC's Cori, it is best to run MODTRAN on the frontend.
To run this tool on the Cori frontend, run the following:
```
make
mpirun -np <procs> mpi_modtran.x input output_dir data_dir
# input is the input JSON file, output_dir is the output directory, and data_dir is the data directory of your MODTRAN installation
```

### FV3_destruct_output.py
This tool is used to destruct MODTRAN's output into a netCDF file.
```
usage: python3 FV3_destruct_output.py [-h] input output cases szas

Destruct FV3 ADM MODTRAN output into netCDF file.

positional arguments:
  input       Path to directory containing MODTRAN's output.
  output      Output file path.
  cases       Number of cases.
  szas        Number of szas.

optional arguments:
  -h, --help  show this help message and exit
```

### print_netcdf.py
This tool is not strictly part of the pipeline. It is a utility to print the structure of netCDF files. It can also be used to print out individual variables.
```
usage: python3 print_netcdf.py [-h] [--var [VAR]] input

Quickly print netCDF file info.

positional arguments:
  input        Input file path.

optional arguments:
  -h, --help   show this help message and exit
  --var [VAR]  Variable to print contents of.
```

### diff_netcdf.py
This tool is not strictly part of the pipeline. It is a utility to print differences in two netCDF files.
```
usage: python3 diff_netcdf.py [-h] input1 input2

Compare 2 netCDF files.

positional arguments:
  input1      First input file path.
  input2      Second input file path.

optional arguments:
  -h, --help  show this help message and exit
```

### graphs/
This directory contains scripts for creating graphs from the output netCDF from the main pipeline. These scripts should be relatively self-explanatory, as the intention is that they are templates for your own creation of plots.
