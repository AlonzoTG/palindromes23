CC=gcc
CFLAGS= -march=native -std=c11 -Wall -lgmp -lgomp

debug:clean 
	$(CC) $(CFLAGS) -ggdb -o palindrome main2.c 

inspect:clean
	$(CC) $(CFLAGS) -S -o palindrome.asm main2.c

stable:clean 
	$(CC) $(CFLAGS) -O3 -ggdb -fopenmp -o palindrome main2.c 

old:clean 
	$(CC) $(CFLAGS) -O3 -fopenmp -o palindrome main.c 

single:clean 
	$(CC) $(CFLAGS) -O3 -ggdb -o palindrome main2.c

clean:
	rm -vfr *~ *.o *.asm palindrome
