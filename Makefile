#MODTRAN_DIR = /nobackup/rarbore/MODTRAN6.0
MODTRAN_DIR = /project/projectdirs/m4098/modtran/modtran/modtran6/mod6_install/MODTRAN6.0

#CXX = mpicxx -std=c++11 -I$(MODTRAN_DIR)/developer/include -L$(MODTRAN_DIR)/bin/linux -lmod6rt -Wall
CXX = g++ -std=c++11 -I$(MODTRAN_DIR)/developer/include -L$(MODTRAN_DIR)/bin/linux -lmod6rt -I$(HOME)/opt/openmpi/build/ompi/include -L/usr/lib64/mpi/gcc/openmpi/lib64 -l:libmpi.so.12 -Wall
RUN = mpirun -quiet

mpi_modtran.x: mpi_modtran.cpp
	$(CXX) -o $@ $^

clean:
	rm -rf *.x *.json *.sap $(SCRATCH)/output/

.PHONY: clean
