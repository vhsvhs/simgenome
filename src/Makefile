CC=g++
MPICC=mpicc
CFLAGS=-Wall
LDFLAGS=
EXECUTABLE=hello

OBJECTS=common.o gene.o fout.o genome.o ga.o fitness.o landscape.o mutate.o population.o psam.o ptable.o serialize.o settings.o main.o

all: simreg

simreg: $(OBJECTS)
	$(CC) $(CFLAGS) $(OBJECTS) -o simreg

simreg-mpi: $(OBJECTS)
	$(MPICC) $(OBJECTS) -o simreg-mpi

ga.o: ga.c
	$(CC) $(CFLAGS) -c ga.c

gene.o: gene.c
	$(CC) $(CFLAGS) -c gene.c

genome.o: genome.c
	$(CC) $(CFLAGS) -c genome.c
	
fout.o: fout.c
	$(CC) $(CFLAGS) -c fout.c

main.o: main.c
	$(CC) $(CFLAGS) -c main.c

scratch.o: scratch.c
	$(CC) $(CFLAGS) -c scratch.c

settings.o: settings.c
	$(CC) $(CFLAGS) -c settings.c
	
serialize.o: serialize.c
	$(CC) $(CFLAGS) -c serialize.c

population.o: population.c
	$(CC) $(CFLAGS) -c population.c

psam.o: psam.c
	$(CC) $(CFLAGS) -c psam.c

ptable.o: ptable.c
	$(CC) $(CFLAGS) -c ptable.c

common.o: common.c
	$(CC) $(CFLAGS) -c common.c

mutate.o: mutate.c
	$(CC) $(CFLAGS) -c mutate.c
	
landscape.o: landscape.c
	$(CC) $(CFLAGS) -c landscape.c

fitness.o: fitness.c
	$(CC) $(CFLAGS) -c fitness.c

clean:
	rm -rf *o hello