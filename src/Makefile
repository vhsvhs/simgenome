CC=gcc
CFLAGS=-c -Wall
LDFLAGS=
SOURCES=common.cpp gene.cpp mutate.cpp psam.cpp scratch.cpp settings.cpp  
OBJECTS=$(SOURCES:.cpp=.o)
EXECUTABLE=hello

all: hello

hello: common.o gene.o fout.o genome.o mutate.o psam.o settings.o scratch.o
	$(CC) common.o gene.o fout.o genome.o mutate.o psam.o settings.o scratch.o -o scratch

gene.o: gene.cpp
	$(CC) -c gene.cpp

genome.o: genome.cpp
	$(CC) -c genome.cpp
	
fout.o: fout.cpp
	$(CC) -c fout.cpp

scratch.o: scratch.cpp
	$(CC) -c scratch.cpp

settings.o: settings.cpp
	$(CC) -c settings.cpp

psam.o: psam.cpp
	$(CC) -c psam.cpp

common.o: common.cpp
	$(CC) -c common.cpp

mutate.o: mutate.cpp
	$(CC) -c mutate.cpp

clean:
	rm -rf *o hello