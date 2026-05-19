ICC=icc
CC=gcc
GCC=gcc

#CC=$(GCC)
CC=$(GCC)
#CC=clang

CHECK=-DCHECK
#SIMPLESTENCIL=-Dsimplestencil

ifeq ($(CC), $(ICC))
	OPT_FLAGS     := -O3 -xHost -ansi-alias -ipo -fp-model precise $(CHECK) $(SIMPLESTENCIL) #-xAVX2 -mfma
	OMP_FLAGS     := -qopenmp
else ifeq ($(CC), $(GCC))
	OPT_FLAGS     := -O3 -march=native -mtune=native -mavx512f $(CHECK) $(SIMPLESTENCIL) #-mavx2 -mfma
	OMP_FLAGS     := -fopenmp
else
	OPT_FLAGS     := -O3 -march=native -mtune=native
	PAR_FLAGS     := 
	OMP_FLAGS     := -fopenmp -I/usr/lib/gcc/x86_64-linux-gnu/13/include
endif