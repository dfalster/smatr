CC=g++
CFLAGS=-O3 -c -Wall
INCS =
LDFLAGS= -lgsl -lgslcblas #For Mac
SOURCES= main.cpp menu.cpp SMAstats.cpp
OBJECTS=$(SOURCES:.cpp=.o)
EXECUTABLE=smatr

all: $(SOURCES) $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS)
	$(CC) $(LIBS) $(LDFLAGS) $(OBJECTS) -o $@

.cpp.o:
	$(CC) $(CFLAGS) $(INCS) $< -o $@

clean:
	rm -f *.o $(EXECUTABLE)
