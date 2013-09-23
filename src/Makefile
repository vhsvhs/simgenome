CC=g++
MPICC=mpicc
CFLAGS=-c -Wall
LDFLAGS=
EXECUTABLE=hello

OBJECTS=common.o gene.o fout.o genome.o ga.o fitness.o landscape.o mutate.o population.o psam.o ptable.o serialize.o settings.o main.o

all: simreg simreg-mpi

simreg: $(OBJECTS)
	$(CC) $(OBJECTS) -o simreg

simreg-mpi: $(OBJECTS)
	$(MPICC) $(OBJECTS) -o simreg-mpi

ga.o: ga.c
	$(CC) -c ga.c

gene.o: gene.c
	$(CC) -c gene.c

genome.o: genome.c
	$(CC) -c genome.c
	
fout.o: fout.c
	$(CC) -c fout.c

main.o: main.c
	$(CC) -c main.c

scratch.o: scratch.c
	$(CC) -c scratch.c

settings.o: settings.c
	$(CC) -c settings.c
	
serialize.o: serialize.c
	$(CC) -c serialize.c

population.o: population.c
	$(CC) -c population.c

psam.o: psam.c
	$(CC) -c psam.c

ptable.o: ptable.c
	$(CC) -c ptable.c

common.o: common.c
	$(CC) -c common.c

mutate.o: mutate.c
	$(CC) -c mutate.c
	
landscape.o: landscape.c
	$(CC) -c landscape.c

fitness.o: fitness.c
	$(CC) -c fitness.c

clean:
	rm -rf *o hello