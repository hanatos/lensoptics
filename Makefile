CC=gcc
CXX=g++
CFLAGS=-std=c11 -Wall -g -Iext/ -D_GNU_SOURCE -fopenmp
CXXFLAGS=-Wall -g -Iext/ -D_GNU_SOURCE -fopenmp
# OPTFLAGS=-I. -ggdb3 -Isrc/
OPTFLAGS=-O3 -ffast-math -mfpmath=sse -march=native -fexpensive-optimizations -DNDEBUG -fno-finite-math-only -I. -Isrc/
LDFLAGS=-lm

HEADERS=\
src/lenssystem.h\
src/poly.h\
src/raytrace_draw.h\
src/raytrace.h\
src/gencode.h\
src/spectrum.h

.PHONY=all clean

all: view fit gencode fresnel simplify

view: Makefile src/view.c ${HEADERS}
	${CC} ${OPTFLAGS} ${CFLAGS} src/view.c $(shell pkg-config --cflags --libs gtk+-2.0) ${LDFLAGS} -o view ${LDFLAGS}

fresnel: Makefile src/fresnel.c ${HEADERS}
	${CC} ${OPTFLAGS} ${CFLAGS} src/fresnel.c -o fresnel ${LDFLAGS}

fit: Makefile src/fit.c ${HEADERS}
	${CXX} ${OPTFLAGS} ${CXXFLAGS} src/fit.c -o fit ${LDFLAGS}

genpoly: Makefile src/genpoly.c ${HEADERS}
	${CC} ${OPTFLAGS} ${CFLAGS} src/genpoly.c -o genpoly ${LDFLAGS}

parsepoly: Makefile src/parsepoly.c ${HEADERS}
	${CC} ${OPTFLAGS} ${CFLAGS} src/parsepoly.c -o parsepoly ${LDFLAGS}

printpoly: Makefile src/printpoly.c ${HEADERS}
	${CC} ${OPTFLAGS} ${CFLAGS} src/printpoly.c -o printpoly ${LDFLAGS}

simplify: Makefile src/simplify.c ${HEADERS}
	${CC} ${OPTFLAGS} ${CFLAGS} src/simplify.c -o simplify ${LDFLAGS}

gencode: Makefile src/gencode.c ${HEADERS}
	${CC} ${OPTFLAGS} ${CFLAGS} src/gencode.c -o gencode ${LDFLAGS}

clean:
	rm -f view fit gencode fresnel fresnel.dat simplify
