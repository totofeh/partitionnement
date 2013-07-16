CC= g++
CFLAGS= -g 
LFLAGS= 
SRCDIR= ./src
SRC= main_gggp.cpp gggp.cpp util.cpp
OBJ= $(SRC:.cpp=.o)

all: gggp

gggp : $(OBJ)
	$(CC) $(CFLAGS) $(OBJ) -o $@ $(LFLAGS)

%.o: $(SRCDIR)/%.cpp 
	$(CC) $(CFLAGS) -c $<

clean:
	rm -f *.o gggp