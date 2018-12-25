CC=g++ -std=c++11
#CC=g++ -std=c++0x

SOURCE   = /Users/jiewang/Google\ Drive/jie_programs
#SOURCE   = /home/jiew/jie_programs

CLIBRARY = $(SOURCE)/JieCLib
EIGEN    = $(SOURCE)/eigen
#SPECTRA  = $(SOURCE)/spectra/include/
SPA_CAL  = $(SOURCE)/Jie_Sparse_Diag/link_to_arpack
DIAG     = $(SOURCE)/FQH_on_Cylinder
CFLAGS   = -I$(CLIBRARY) -I$(EIGEN) -I$(DIAG) 

a.out:  main.o quantumscar.o
	$(CC) -O3 $(CFLAGS) -fopenmp -o result main.o quantumscar.o $(DIAG)/pseu_pot_diag.o $(DIAG)/ModelWf.o $(DIAG)/phpfaffian.o $(CLIBRARY)/utility.o $(CLIBRARY)/weir3.o $(CLIBRARY)/z_function_m.o -lgfortran

clean:
	rm -f *~ *.o a.out

%.o:	%.cpp
	$(CC) -w -O3 -fopenmp $(CFLAGS) -c $<

%.o:	%.f90
	gfortran -w -O3 $(CFLAGS) -c $<
