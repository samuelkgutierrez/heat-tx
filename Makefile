all: heat-tx

CFLAGS = -Wall -O3

heat-tx: heat-tx.c

clean:
	rm -f heat-tx
	rm -rf heat-tx.dSYM
