#MODTRAN_DIR = /nobackup/rarbore/MODTRAN6.0
MODTRAN_DIR = /project/projectdirs/m4098/modtran/modtran/modtran6/mod6_install/MODTRAN6.0
OPT_DIR = /project/projectdirs/m4098/libRadtran-install
LIBRADTRAN_INSTALL = /project/projectdirs/m4098/libRadtran-install/libRadtran/bin
OPENMPI_INSTALL = /project/projectdirs/m4098/libRadtran-install/openmpi

#CXX = mpicxx -std=c++11 -I$(MODTRAN_DIR)/developer/include -L$(MODTRAN_DIR)/bin/linux -lmod6rt -Wall
CXX = g++ -std=c++11 -I$(MODTRAN_DIR)/developer/include -L$(MODTRAN_DIR)/bin/linux -lmod6rt -I$(OPENMPI_INSTALL)/build/ompi/include -L/usr/lib64/mpi/gcc/openmpi/lib64 -l:libmpi.so.12 -Wall
RUN = mpirun -quiet

all: mpi_modtran.x wc.gamma_007.0.mie.cdf

mpi_modtran.x: mpi_modtran.cpp
	$(CXX) -o $@ $^
wc.gamma_007.0.mie.cdf: mie_input.txt
	LD_LIBRARY_PATH=$(LD_LIBRARY_PATH):$(OPT_DIR)/libRadtran/lib:$(OPT_DIR)/gsl/lib:$(OPT_DIR)/netcdf/lib:$(OPT_DIR)/hdf5/lib $(LIBRADTRAN_INSTALL)/mie <mie_input.txt

clean:
	rm -rf *.x *.json *.sap *.cdf $(SCRATCH)/output/

.DEFAULT: all

