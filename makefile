#CFLAGS=-Wall -g -I/usr/local/opt/readline/include `pkg-config --cflags gsl` -std=c99
CFLAGS=-Wall -g `pkg-config --cflags gsl`  -std=c99
#LDFLAGS=-g -L/usr/local/opt/readline/lib `pkg-config --libs gsl` -lm 
LDFLAGS=-g `pkg-config --libs gsl` -lreadline -lm 

current_dir = $(shell pwd)

all: libodexp.dylib

libodexp.dylib: odexp.o methods_odexp.o
		gcc $(LDFLAGS) -dynamiclib -install_name $(current_dir)/libodexp.dylib   -lm -o libodexp.dylib odexp.o methods_odexp.o 

odexp.o: odexp.c odexp.h methods_odexp.h
	    gcc -c odexp.c $(CFLAGS) 

methods_odexp.o: methods_odexp.c methods_odexp.h odexp.h
		gcc -c methods_odexp.c $(CFLAGS)


clean:
		rm -f *.o

veryclean:
		rm -f *.o; rm -f *.out; rm -rf *.out.dSYM; rm *.dylib
