CFLAGS = -fopenmp -O3 -std=c99

all: symmetrize xpose

symmetrize: clean
	gcc $(CFLAGS) symmetrize.c -o ./bin/symmetrize

xpose: clean
	gcc $(CFLAGS) xpose.c -o ./bin/xpose


clean:
	rm -f ./bin/*
