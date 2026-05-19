EXE=exe-3d-jacobi-27p

include ../common.make

all: $(EXE)

$(EXE): define.h	main.c		naive_scalar.c	naive_vector.c	vectime_transpose_boundary_extra_array.c	vectime_transpose_boundary_extra_array_512.c
	$(CC) $(OPT_FLAGS) $^ -lm -o $(EXE)

clean:
	rm -f  exe*
