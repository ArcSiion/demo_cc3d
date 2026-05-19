EXE=exe-3d-cc3d-diffusion

include ../common.make

all: $(EXE)

$(EXE): *.c *.h
	$(CC) $(OPT_FLAGS) $(OMP_FLAGS) $^ -o $(EXE) -lm

clean:
	rm -f exe*
