ncincdir = -I/home/public/easybuild/software/netCDF-Fortran/4.5.3-gmpich-2021.01/include
nclibdir = -L/home/public/easybuild/software/netCDF-Fortran/4.5.3-gmpich-2021.01/lib
nclib    = -lnetcdff

mods = overprint.o         \
       common_vars.o       \
       netcdf_error.o      \
       esatdesdT.o         \
       coordsmod.o         \
       read_data.o         \
       airmass.o           \
       daily.o             \
       daily_insolation.o  \
       netrad_pet.o        \
       orbit.o             \
       crop_suit.o         \
       surfrad.o           \
       aet_alpha.o         \
       rad_and_pet.o       \
       biome1mod.o

biome1 = biome1driver.o
pointdriver = biome1pointdriver.o

FC = gfortran
CC = gcc
FCFLAGS = -ffree-line-length-none $(ncincdir)
CFLAGS =

%.mod : %.o
	$(MAKE) $<

%.o : %.f90
	$(FC) $(FCFLAGS) -c -o $@ $<

biome1 ::	$(mods) $(biome1)
	$(FC) $(FCFLAGS) -o biome1 $(biome1) $(mods) $(nclibdir) $(nclib)

point ::	$(mods) $(pointdriver)
	$(FC) $(FCFLAGS) -o biome1point $(pointdriver) $(mods) $(nclibdir) $(nclib)

all :: biome1

clean :: 
	rm *.mod *.o biome1
