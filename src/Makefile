CC=gcc
CFLAGS=-c -Wall
LDFLAGS=
SOURCES=common.cpp mutate.cpp psam.cpp scratch.cpp settings.cpp  
OBJECTS=$(SOURCES:.cpp=.o)
EXECUTABLE=hello

all: hello

hello: common.o mutate.o psam.o settings.o scratch.o
	$(CC) common.o mutate.o psam.o settings.o scratch.o -o scratch

mutate.o: mutate.cpp
	$(CC) -c mutate.cpp

scratch.o: scratch.cpp
	$(CC) -c scratch.cpp

settings.o: settings.cpp
	$(CC) -c settings.cpp

psam.o: psam.cpp
	$(CC) -c psam.cpp

common.o: common.cpp
	$(CC) -c common.cpp

clean:
	rm -rf *o hello