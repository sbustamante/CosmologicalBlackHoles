CC=gcc
CFLAGS = -Wall -I/home/bustamsn/Software/gsl-2.1/include
LDFLAGS = -L/home/bustamsn/Software/gsl-2.1/lib
LIBS= -lm -lgsl -lgslcblas
SOURCES=main.c
OBJECTS= main.o
EXECUTABLE=main

$(EXECUTABLE): $(OBJECTS)
	$(CC) $(LDFLAGS) $(OBJECTS) $(LIBS) -o $(EXECUTABLE)
	
clean:
	rm -f $(OBJECTS) $(EXECUTABLE)
