#fort = mpif90 -check all
fort = mpiifort -lfftw3 -lfftw3_threads 
RM = rm -f
OBJ =  main.o  funct.o makegrid.o ksum.o conv.o fftw3.o

../MO_IPT.run : $(OBJ)
	$(fort) -O3 -o ../MO_IPT.run $(OBJ) 

main.o : main.f90
	$(fort) -O3 -c -o main.o main.f90

makegrid.o : makegrid.f90 
	$(fort) -O3 -c -o makegrid.o makegrid.f90 

funct.o : funct.f 
	$(fort) -O3 -c -o funct.o funct.f 

ksum.o : ksum.f90
	$(fort) -O3 -c -o ksum.o ksum.f90

conv.o : conv.f90
	$(fort) -O3 -c -o conv.o conv.f90

fftw3.o : fftw3.f90
	$(fort) -O2 -c -o fftw3.o fftw3.f90

clean:
	$(RM) $(OBJ) *mod ../MO_IPT.run 
