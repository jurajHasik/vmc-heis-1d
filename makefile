FC=gfortran

MT_FLAGS=-fno-range-check

MT=prng/MersenneTwister.f90
# Suppres warning during compilation of MersenneTwister.f90
# See https://gcc.gnu.org/onlinedocs/gfortran/Fortran-Dialect-Options.html
MT_FLAGS=-fno-range-check

vmc: clean
	$(FC) $(MT_FLAGS) -O3 -o vmc.x $(MT) vmc-heis-1d.f90
	rm -f *.mod
clean:
	rm -f *.mod *.x