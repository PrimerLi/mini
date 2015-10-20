cc = gfortran
lib = -L/usr/local/lib -llapack -lblas
obj = Chi.o inverse.o 
all: Chi.out
Chi.o:\
        $(cc) -c Chi.f90
inverse.o:
	$(cc) -c inverse.f90
Chi.out: $(obj)
	$(cc) -o Chi.out $(obj) $(lib)
clean:
	rm Chi.out $(obj)
