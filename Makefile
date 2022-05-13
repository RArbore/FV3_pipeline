#MODTRAN_DIR = /nobackup/rarbore/MODTRAN6.0
MODTRAN_DIR = /project/projectdirs/m4098/modtran/modtran/modtran6/mod6_install/MODTRAN6.0

#CXX = mpicxx -std=c++11 -I$(MODTRAN_DIR)/developer/include -L$(MODTRAN_DIR)/bin/linux -lmod6rt -Wall
CXX = g++ -std=c++11 -I$(MODTRAN_DIR)/developer/include -L$(MODTRAN_DIR)/bin/linux -lmod6rt -I$(HOME)/opt/openmpi/build/ompi/include -L/usr/lib64/mpi/gcc/openmpi/lib64 -l:libmpi.so.12 -Wall
RUN = mpirun -quiet

all: mpi_modtran.x wc.gamma_007.0.mie.cdf

mpi_modtran.x: mpi_modtran.cpp
	$(CXX) -o $@ $^
wc.gamma_007.0.mie.cdf: mie_input.txt
	mie <mie_input.txt

clean:
	rm -rf *.x *.json *.sap $(SCRATCH)/output/

.DEFAULT: all
.PHONY: clean, all
