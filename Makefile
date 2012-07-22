all: heat-tx

CFLAGS = -Wall -g -O0

heat-tx: heat-tx.c

clean:
	rm -f heat-tx
	rm -rf heat-tx.dSYM
