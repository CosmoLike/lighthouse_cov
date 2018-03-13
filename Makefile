
home:
	gcc  multi_covariances_real_mpp_home.c -o ./multi_covariances_real_mpp_home  -lgsl -lfftw3 -lgslcblas -lm -O3 -std=gnu99 -ffast-math -funroll-loops -L../cosmolike_light/class -lclass

class:
	cd ../cosmolike_light/class; $(MAKE)

