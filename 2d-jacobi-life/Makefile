

EXE=exe-2d-jacobi-life

include ../common.make

all: $(EXE)

$(EXE): main.c naive_scalar.c naive_vector.c vectime.c define.h vectime_extra_array.c vectime_512.c vectime_extra_array_512.c
	$(CC) $(OPT_FLAGS)  $^ -o $(EXE)

clean:
	rm -f  exe* 


