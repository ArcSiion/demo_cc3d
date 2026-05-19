exe = \
	2d-jacobi-life\
	3d-jacobi-7p\
	3d-jacobi-27p\
	3d-cc3d-diffusion

all:
	@-for d in $(exe); do \
		make -C $$d; \
	done

clean:
	@-for d in $(exe); do \
		make -C $$d clean; \
	done
