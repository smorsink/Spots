# Source (C) Coire Cadeau 2007, all rights reserved.

# Permission is granted for private use only, and not
# distribution, either verbatim or of derivative works,
# in whole or in part.

# The code is not thoroughly tested or guaranteed for
# any particular use.

#Sharon needs to use the following:
CC=c++
CCFLAGS=-Wall -pedantic -O3 

LDFLAGS=-lm

NAMES=spot emit bend detect

OBJ=PolyOblModelBase.o  PolyOblModelCFLQS.o PolyOblModelNHQS.o Units.o OblDeflectionTOA.o \
	Chi.o Hydrogen.o Atmo.o TimeDelays.o BlackBody.o Instru.o Ism.o SphericalOblModel.o matpack.o interp.o nrutil.o # defining the objects

EOBJ=PolyOblModelBase.o  PolyOblModelCFLQS.o PolyOblModelNHQS.o Units.o OblDeflectionTOA.o \
	Chi.o Hydrogen.o Atmo.o TimeDelays.o BlackBody.o SphericalOblModel.o matpack.o interp.o nrutil.o # defining the objects

DOBJ=Units.o \
	Instru.o Ism.o matpack.o interp.o nrutil.o # defining the objects

APPOBJ=Spot.o

BOBJ=Bend.o

all: $(NAMES)

spot: Spot.o $(OBJ)
	$(CC) $(CCFLAGS) Spot.o $(OBJ) $(LDFLAGS) -o spot

emit: SpotEmit.o $(EOBJ)
	$(CC) $(CCFLAGS) SpotEmit.o $(EOBJ) $(LDFLAGS) -o emit

detect: SpotDetect.o $(DOBJ)
	$(CC) $(CCFLAGS) SpotDetect.o $(DOBJ) $(LDFLAGS) -o detect

bend: Bend.o $(OBJ)
	$(CC) $(CCFLAGS) Bend.o $(OBJ) $(LDFLAGS) -o bend

Spot.o: \
	Spot.cpp \
	OblDeflectionTOA.h \
	Chi.h \
	Atmo.h \
	TimeDelays.h \
	Instru.h \
	Ism.h \
	Struct.h \
	PolyOblModelNHQS.h \
	PolyOblModelCFLQS.h \
	SphericalOblModel.h \
	OblModelBase.h \
	Units.h \
	Makefile
	$(CC) $(CCFLAGS) -c Spot.cpp

SpotEmit.o: \
	SpotEmit.cpp \
	OblDeflectionTOA.h \
	Chi.h \
	Atmo.h \
	TimeDelays.h \
	Struct.h \
	PolyOblModelNHQS.h \
	PolyOblModelCFLQS.h \
	SphericalOblModel.h \
	OblModelBase.h \
	Units.h \
	Makefile
	$(CC) $(CCFLAGS) -c SpotEmit.cpp

SpotDetect.o: \
	SpotDetect.cpp \
	Struct.h \
	Units.h \
	Makefile
	$(CC) $(CCFLAGS) -c SpotDetect.cpp


Bend.o: \
	Bend.cpp \
	OblDeflectionTOA.h \
	Chi.h \
	Atmo.h \
	Instru.h \
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
	McPhac.h \
	Units.h \
	matpack.h
	$(CC) $(CCFLAGS) -c Atmo.cpp

Hydrogen.o: \
	Hydrogen.h \
	OblDeflectionTOA.h \
	Hydrogen.cpp \
	OblModelBase.h \
	McPhac.h \
	Units.h \
	matpack.h
	$(CC) $(CCFLAGS) -c Hydrogen.cpp

TimeDelays.o: \
	TimeDelays.h \
	TimeDelays.cpp \
	interp.h \
	Units.h \
	matpack.h
	$(CC) $(CCFLAGS) -c TimeDelays.cpp

McPhac.o: \
	McPhac.h \
	McPhac.cpp \
	interp.h \
	Units.h \
	matpack.h
	$(CC) $(CCFLAGS) -c McPhac.cpp

BlackBody.o: \
	BlackBody.h \
	BlackBody.cpp \
	interp.h \
	Units.h \
	matpack.h
	$(CC) $(CCFLAGS) -c BlackBody.cpp

Instru.o: \
	Instru.h \
	OblDeflectionTOA.h \
	Instru.cpp \
	OblModelBase.h \
	Units.h \
	matpack.h
	$(CC) $(CCFLAGS) -c Instru.cpp

Ism2.o: \
	Ism.h \
	Ism2.cpp \
	Units.h 
	$(CC) $(CCFLAGS) -c Ism2.cpp


Instru2.o: \
	Instru.h \
	OblDeflectionTOA.h \
	Instru2.cpp \
	OblModelBase.h \
	Units.h \
	matpack.h
	$(CC) $(CCFLAGS) -c Instru2.cpp

Ism.o: \
	Ism.h \
	Ism.cpp \
	Units.h 
	$(CC) $(CCFLAGS) -c Ism.cpp



Units.o: \
	Units.h \
	Units.cpp
	$(CC) $(CCFLAGS) -c Units.cpp

interp.o: \
	interp.h \
	nrutil.h \
	interp.cpp
	$(CC) $(CCFLAGS) -c interp.cpp

nrutil.o: \
	nrutil.h \
	nrutil.c
	$(CC) $(CCFLAGS) -c nrutil.c

matpack.o: \
	matpack.h \
	matpack.cpp \
	Exception.h
	$(CC) $(CCFLAGS) -c matpack.cpp

clean:
	rm -f core *~ $(OBJ) $(APPOBJ)

veryclean:
	rm -f core *~ $(OBJ) $(APPOBJ) $(NAMES)
