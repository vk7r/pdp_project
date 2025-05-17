CC = mpicc
CFLAGS = -O3 #-Wall
BINS = main
# -fsanitize=address -g

all: $(BINS)

main: matrix.h matrix.c power_method.h power_method.c main.c
	$(CC) $(CFLAGS) -o power_method power_method.c matrix.c main.c -lm

clean:
	$(RM) power_method
