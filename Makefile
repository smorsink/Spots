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

NAMES=spot bend

OBJ=PolyOblModelBase.o  PolyOblModelCFLQS.o PolyOblModelNHQS.o Units.o OblDeflectionTOA.o \
	Chi.o Atmo.o SphericalOblModel.o matpack.o # defining the objects

APPOBJ=Spot.o

BOBJ=Bend.o

all: $(NAMES)

spot: Spot.o $(OBJ)
	$(CC) $(CCFLAGS) Spot.o $(OBJ) $(LDFLAGS) -o spot

bend: Bend.o $(OBJ)
	$(CC) $(CCFLAGS) Bend.o $(OBJ) $(LDFLAGS) -o bend

Spot.o: \
	Spot.cpp \
	OblDeflectionTOA.h \
	Chi.h \
	Atmo.h \
	Struct.h \
	PolyOblModelNHQS.h \
	PolyOblModelCFLQS.h \
	SphericalOblModel.h \
	OblModelBase.h \
	Units.h \
	Makefile
	$(CC) $(CCFLAGS) -c Spot.cpp

Bend.o: \
	Bend.cpp \
	OblDeflectionTOA.h \
	Chi.h \
	Atmo.h \
	Struct.h \
	PolyOblModelNHQS.h \
	PolyOblModelCFLQS.h \
	SphericalOblModel.h \
	OblModelBase.h \
	Units.h \
	Makefile
	$(CC) $(CCFLAGS) -c Bend.cpp

PolyOblModelBase.o: \
	PolyOblModelBase.h \
	PolyOblModelBase.cpp \
	OblModelBase.h
	$(CC) $(CCFLAGS) -c PolyOblModelBase.cpp

PolyOblModelCFLQS.o: \
	PolyOblModelCFLQS.h \
	PolyOblModelCFLQS.cpp \
	PolyOblModelBase.h
	$(CC) $(CCFLAGS) -c PolyOblModelCFLQS.cpp

PolyOblModelNHQS.o: \
	PolyOblModelNHQS.h \
	PolyOblModelNHQS.cpp \
	PolyOblModelBase.h
	$(CC) $(CCFLAGS) -c PolyOblModelNHQS.cpp


SphericalOblModel.o: \
	SphericalOblModel.h \
	SphericalOblModel.cpp \
	OblModelBase.h
	$(CC) $(CCFLAGS) -c SphericalOblModel.cpp

OblDeflectionTOA.o: \
	OblDeflectionTOA.h \
	OblDeflectionTOA.cpp \
	OblModelBase.h \
	Units.h \
	matpack.h
	$(CC) $(CCFLAGS) -c OblDeflectionTOA.cpp


Chi.o: \
	Chi.h \
	OblDeflectionTOA.h \
	Chi.cpp \
	OblModelBase.h \
	Units.h \
	matpack.h
	$(CC) $(CCFLAGS) -c Chi.cpp

Atmo.o: \
	Atmo.h \
	OblDeflectionTOA.h \
	Atmo.cpp \
	OblModelBase.h \
	Units.h \
	matpack.h
	$(CC) $(CCFLAGS) -c Atmo.cpp


Units.o: \
	Units.h \
	Units.cpp
	$(CC) $(CCFLAGS) -c Units.cpp

matpack.o: \
	matpack.h \
	matpack.cpp \
	Exception.h
	$(CC) $(CCFLAGS) -c matpack.cpp

clean:
	rm -f core *~ $(OBJ) $(APPOBJ)

veryclean:
	rm -f core *~ $(OBJ) $(APPOBJ) $(NAMES)
