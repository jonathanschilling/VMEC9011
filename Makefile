
all: vmec

vmec: crlib64.f code.f vac.f name0 name1 name2 name3 name4
	gfortran vac.f -c -o vac.f.o
	gfortran -std=legacy crlib64.f -c -o crlib64.f.o
	gfortran -fbacktrace -g -fdefault-real-8 -std=legacy code.f -c -o code.f.o
	gfortran code.f.o crlib64.f.o vac.f.o -o vmec
