# Source (C) Coire Cadeau 2007, all rights reserved.

# Permission is granted for private use only, and not
# distribution, either verbatim or of derivative works,
# in whole or in part.

# The code is not thoroughly tested or guaranteed for
# any particular use.

#Albert's version of the flags
#CC=g++
#CCFLAGS=-Wall -pedantic -O3 -std=c++11

#Sharon needs to use the following:
CC=c++
CCFLAGS=-Wall -pedantic -O3 

LDFLAGS=-lm

NAMES=format

OBJ=Format.o

all: $(NAMES)

format: Format.o
	$(CC) $(CCFLAGS) Format.o $(LDFLAGS) -o format

Format.o: \
	Format.cpp \
	Makefile
	$(CC) $(CCFLAGS) -c Format.cpp


clean:
	rm -f core *~ $(OBJ) $(APPOBJ)

veryclean:
	rm -f core *~ $(OBJ) $(APPOBJ) $(NAMES)
