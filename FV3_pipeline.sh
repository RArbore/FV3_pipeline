#!/usr/bin/sh
set -e

if [ $# -eq 0 ]; then
	>&2 echo "Please provide an input netCDF file with the expected data fields."
	exit
fi

echo "Please ignore all errors after this point pertaining to invalid arguments. Scripts called after this point have their own semantics for arguments, and are not indicative of the arguments you should use for this script. If you receive argument errors for any scripts after this point, it is likely a bug in this script."

DATA_DIR=/project/projectdirs/m4098/modtran/modtran/modtran6/mod6_install/MODTRAN6.0/DATA
BIN_DIR=/project/projectdirs/m4098/modtran/modtran/modtran6/mod6_install/MODTRAN6.0/bin/linux
OPT_DIR=/project/projectdirs/m4098/libRadtran-install

MIE_INPUT='mie_input.txt'
MIE_OUTPUT='wc.gamma_007.0.mie.cdf'
JSON_FILE=$(echo $1 | sed 's/.nc/.json/')
MODTRAN_FILE=$(echo $1 | sed 's/.nc/_modtran_input.json/')
SAP_FILE=$(echo $1 | sed 's/.nc/.sap/')

OUTPUT_DIR=$SCRATCH/output/
OUTPUT_FILE=$(echo $1 | sed 's/.nc/_modtran_output.nc/')
LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$BIN_DIR:/usr/lib64/mpi/gcc/openmpi/lib64:$OPT_DIR/libRadtran/lib:$OPT_DIR/gsl/lib:$OPT_DIR/netcdf/lib:$OPT_DIR/hdf5/lib

make all
python3 json_netcdf.py $1 $JSON_FILE
python3 FV3_construct_input.py $JSON_FILE $SAP_FILE $MODTRAN_FILE
python3 generate_sap.py $MIE_OUTPUT $MODTRAN_FILE $SAP_FILE
mkdir -p $OUTPUT_DIR
mpirun -np 32 mpi_modtran.x $MODTRAN_FILE $OUTPUT_DIR $DATA_DIR
python3 FV3_destruct_output.py $OUTPUT_DIR $OUTPUT_FILE 90 10
