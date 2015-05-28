CFLAGS=-Wall -g `pkg-config --cflags gsl` -pedantic
LDFLAGS=-g `pkg-config --libs gsl` -lm

all: libodexp.dylib

libodexp.dylib: odexp.o
		gcc $(LDFLAGS)  -dynamiclib -install_name ~/Documents/codes/odexp/libodexp.dylib -lm -o libodexp.dylib odexp.o 

odexp.o: odexp.c odexp.h
	    gcc -c odexp.c $(CFLAGS) 


clean:
		rm -f *.o

veryclean:
		rm -f *.o; rm -f *.out; rm -rf *.out.dSYM
