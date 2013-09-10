CC=gcc
CFLAGS=-c -Wall
LDFLAGS=
SOURCES=scratch.cpp psam.cpp common.cpp
OBJECTS=$(SOURCES:.cpp=.o)
EXECUTABLE=hello

all: hello

hello: common.o psam.o scratch.o
	$(CC) common.o psam.o scratch.o -o scratch

scratch.o: scratch.cpp
	$(CC) -c scratch.cpp

psam.o: psam.cpp
	$(CC) -c psam.cpp

common.o: common.cpp
	$(CC) -c common.cpp

clean:
	rm -rf *o hello