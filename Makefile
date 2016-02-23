CC=gcc
CFLAGS=-std=c11 -Wall -g -I../ext/levmar-2.6 -D_GNU_SOURCE
OPTFLAGS=-I. -ggdb3
#OPTFLAGS=-O3 -ffast-math -mfpmath=sse -march=native -fexpensive-optimizations -DNDEBUG -fno-finite-math-only
LD_LEVMAR=../ext/levmar-2.6/liblevmar.a -llapack -lblas
LDFLAGS=-lm

HEADERS=\
lenssystem.h\
poly.h\
raytrace_draw.h\
raytrace.h\
spectrum.h

.PHONY=all clean

all: view fit

view: Makefile view.c ${HEADERS}
	${CC} ${OPTFLAGS} ${CFLAGS} view.c $(shell pkg-config --cflags --libs gtk+-2.0) ${LDFLAGS} -o view ${LDFLAGS}

fit: Makefile fit.c ${HEADERS}
	${CC} ${OPTFLAGS} ${CFLAGS} fit.c -o fit ${LDFLAGS} ${LD_LEVMAR}

clean:
	rm -f view fit
