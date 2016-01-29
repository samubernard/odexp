CFLAGS=-Wall -g `pkg-config --cflags gsl`  -pedantic
LDFLAGS=-g `pkg-config --libs gsl` -lreadline -lm 

all: libodexp.dylib

libodexp.dylib: odexp.o methods_odexp.o
		gcc $(LDFLAGS) -dynamiclib -install_name $(HOME)/Documents/codes/odexp/libodexp.dylib -lm -o libodexp.dylib odexp.o methods_odexp.o 

odexp.o: odexp.c odexp.h methods_odexp.h
	    gcc -c odexp.c $(CFLAGS) 

methods_odexp.o: methods_odexp.c methods_odexp.h odexp.h
		gcc -c methods_odexp.c $(CFLAGS)


clean:
		rm -f *.o

veryclean:
		rm -f *.o; rm -f *.out; rm -rf *.out.dSYM; rm *.dylib
