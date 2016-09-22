CC = gcc
CFLAGS = -fopenmp -O3 -std=c99
LDFLAGS = -llapack -lblas -lm

all: cosine symmetrize xpose

cosine: clean
	$(CC) $(CFLAGS) ../dense.c cosine.c $(LDFLAGS) -o ./bin/cosine

symmetrize: clean
	$(CC) $(CFLAGS) symmetrize.c -o ./bin/symmetrize

xpose: clean
	$(CC) $(CFLAGS) xpose.c -o ./bin/xpose


clean:
	rm -f ./bin/*
