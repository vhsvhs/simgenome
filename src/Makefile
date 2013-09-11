CC=gcc
CFLAGS=-c -Wall
LDFLAGS=
SOURCES=common.cpp gene.cpp mutate.cpp psam.cpp scratch.cpp settings.cpp  
OBJECTS=$(SOURCES:.cpp=.o)
EXECUTABLE=hello

all: scratch

scratch: common.o gene.o fout.o genome.o fitness.o landscape.o mutate.o population.o psam.o ptable.o settings.o scratch.o
	$(CC) common.o gene.o fout.o genome.o fitness.o landscape.o mutate.o population.o psam.o ptable.o settings.o scratch.o -o scratch

gene.o: gene.cpp
	$(CC) -c gene.cpp

genome.o: genome.cpp
	$(CC) -c genome.cpp
	
fout.o: fout.cpp
	$(CC) -c fout.cpp

main.o: main.cpp
	$(CC) -c main.cpp

scratch.o: scratch.cpp
	$(CC) -c scratch.cpp

settings.o: settings.cpp
	$(CC) -c settings.cpp

population.o: population.cpp
	$(CC) -c population.cpp

psam.o: psam.cpp
	$(CC) -c psam.cpp

ptable.o: ptable.cpp
	$(CC) -c ptable.cpp

common.o: common.cpp
	$(CC) -c common.cpp

mutate.o: mutate.cpp
	$(CC) -c mutate.cpp
	
landscape.o: landscape.cpp
	$(CC) -c landscape.cpp

fitness.o: fitness.cpp
	$(CC) -c fitness.cpp

clean:
	rm -rf *o hello