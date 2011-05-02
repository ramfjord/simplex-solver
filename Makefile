CC = gcc
#CFLAGS = -std=gnu99 `pkg-config --libs sndfile` `pkg-config  --libs fftw3`
CFLAGS = -std=c99 -lm -lreadline 

simplex: simplex.c matrix.o
	$(CC) $(CFLAGS) -o simplex simplex.c matrix.o

clean:
	rm -rf simplex *.o
