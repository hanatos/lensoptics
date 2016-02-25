CC=gcc
CFLAGS=-std=c11 -Wall -g -Iext/levmar-2.6 -D_GNU_SOURCE -fopenmp
# OPTFLAGS=-I. -ggdb3 -Isrc/
OPTFLAGS=-O3 -ffast-math -mfpmath=sse -march=native -fexpensive-optimizations -DNDEBUG -fno-finite-math-only -I. -Isrc/
LD_LEVMAR=ext/levmar-2.6/liblevmar.a -llapack -lblas
LDFLAGS=-lm

HEADERS=\
src/lenssystem.h\
src/poly.h\
src/raytrace_draw.h\
src/raytrace.h\
src/gencode.h\
src/spectrum.h

.PHONY=all clean

all: view fit gencode fresnel

view: Makefile src/view.c ${HEADERS}
	${CC} ${OPTFLAGS} ${CFLAGS} src/view.c $(shell pkg-config --cflags --libs gtk+-2.0) ${LDFLAGS} -o view ${LDFLAGS}

fresnel: Makefile src/fresnel.c ${HEADERS}
	${CC} ${OPTFLAGS} ${CFLAGS} src/fresnel.c -o fresnel ${LDFLAGS}

ext/levmar-2.6/liblevmar.a:
	make -C ext/levmar-2.6 liblevmar.a

fit: Makefile src/fit.c ${HEADERS} ext/levmar-2.6/liblevmar.a
	${CC} ${OPTFLAGS} ${CFLAGS} src/fit.c -o fit ${LDFLAGS} ${LD_LEVMAR}

gencode: Makefile src/gencode.c ${HEADERS}
	${CC} ${OPTFLAGS} ${CFLAGS} src/gencode.c -o gencode ${LDFLAGS}

clean:
	rm -f view fit gencode fresnel fresnel.dat
