CFLAGS=-Wall -I/usr/local/include -std=c99
#CFLAGS=-Wall -g `pkg-config --cflags gsl`  -std=c99
LDFLAGS=-L/usr/local/lib -lreadline -lgsl -lgslcblas -lm 
#LDFLAGS=-g `pkg-config --libs gsl` -lreadline -lm 

current_dir = $(shell pwd)

all: libodexp.dylib

libodexp.dylib: odexp.o methods_odexp.o  utils_odexp.o rand_gen.o
		gcc $(LDFLAGS) -dynamiclib -install_name $(current_dir)/libodexp.dylib   -lm -o libodexp.dylib odexp.o methods_odexp.o  utils_odexp.o rand_gen.o

odexp.o: odexp.c odexp.h methods_odexp.h
	    gcc -c odexp.c $(CFLAGS) 

methods_odexp.o: methods_odexp.c methods_odexp.h odexp.h
		gcc -c methods_odexp.c $(CFLAGS)

utils_odexp.o: utils_odexp.c utils_odexp.h odexp.h
		gcc -c utils_odexp.c $(CFLAGS)

rand_gen.o: rand_gen.c rand_gen.h odexp.h
		gcc -c rand_gen.c $(CFLAGS)

clean:
		rm -f *.o

veryclean:
		rm -f *.o; rm -f *.out; rm -rf *.out.dSYM; rm *.dylib
