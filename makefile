

NICE    = ifort
Gfort   = gfortran44
FCOMP   = mpif90
F77     = gfortran
CCOMP   = mpicc
CXXCOMP = mpicxx

HYPRE_DIR = /usr/local/lib #HYPRE_DIR = /afs/crc.nd.edu/x86_64_linux/scilib/hypre/2.0.0/intel
HYPRE_LIBS =  -L$(HYPRE_DIR)-lHYPRE #-lHYPRE_struct_ls -lHYPRE_struct_mv -lHYPRE_krylov -lHYPRE_utilities -lHYPRE_DistributedMatrix -lHYPRE_IJ_mv

FOPTS   = -shared-intel -mcmodel=medium
CFLAGS  = -O3
#-----------------------------------------------------------------------------------------------

EXE = Mon_Contact_May26_75.4

OBJ =  M_Math.o     M_General.o     M_Mesh.o   M_Wave_Gen.o     M_SolidFinder_Num.o   M_SolidFinder_Anal.o     M_Solid.o     M_Tower.o       M_Tether.o    M_Add_Force.o     Code_101_O.o 
INC =  M_Math.f90   M_General.f90   M_Mesh.f90 M_Wave_Gen.f90   M_SolidFinder_Num.f90 M_SolidFinder_Anal.f90   M_Solid.f90   M_Tower.f90     M_Tether.f90  M_Add_Force.f90   Code_101_O.f90 

$(EXE): $(OBJ)
	$(NICE) $(FOPTS) -o $(EXE) $(OBJ) $(CFLAGS)

$(OBJ): $(INC)
        
%.o:	%.f90 
	$(NICE) $(FOPTS) -c $<  $(CFLAGS)

clean:
	rm -f *.o *.mod *__genmod.f90

