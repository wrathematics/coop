all: symmetrize xpose

symmetrize: clean
	gcc -fopenmp -O3 symmetrize.c -o symmetrize

xpose: clean
	gcc -fopenmp -O3 xpose.c -o xpose


clean:
	rm -f symmetrize xpose
