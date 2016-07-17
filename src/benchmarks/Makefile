CFLAGS = -fopenmp -O3 -std=c99

all: symmetrize xpose

symmetrize: clean
	gcc $(CFLAGS) symmetrize.c -o symmetrize

xpose: clean
	gcc $(CFLAGS) xpose.c -o xpose


clean:
	rm -f symmetrize xpose
