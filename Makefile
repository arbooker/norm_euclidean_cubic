CC     = gcc
CFLAGS = -O2 -march=native

all: cubic

cubic: cubic.c
	$(CC) $(CFLAGS) -o $@ $< -lm

clean:
	rm -f cubic

.PHONY: all clean
