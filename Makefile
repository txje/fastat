CC=gcc
CFLAGS=-O2
#-std=c99 

OBJECTS = fastat

all: $(OBJECTS)

fastat: fastat.c
	$(CC) $(CFLAGS) fastat.c -o fastat -lz

install: fastat
	cp fastat ~/local/bin/
	chmod 744 ~/local/bin/fastat

.PHONY: clean
clean:
	-rm $(OBJECTS)
