# makefile for LGPL-licensed LAStools
#
#COPTS    = -g -Wall -Wno-deprecated -DDEBUG 

COPTS     = -O3 -Wall -Wno-deprecated -Wno-write-strings -Wno-strict-aliasing -Wno-unused-but-set-variable -Wno-unused-result -DNDEBUG 
#COMPILER  = CC
COMPILER ?= mpic++
LINKER ?= mpic++
#BITS     = -64

#LIBS     = -L/usr/lib64
#LIBS     = -L/usr/lib32
#INCLUDE  = -I/usr/include
INCLUDE = 

LASLIBS     = -L../LASlib/lib
LASINCLUDE  = -I../LASlib/zip/src -I../LASlib/inc 

#all: laszip lasinfo lasprecision las2txt txt2las las2las lasdiff lasmerge lasindex
all: plas lasinfo

plas: plas.o geoprojectionconverter.o ../LASlib/lib/liblas.a
	${LINKER} ${BITS} ${COPTS} plas.o geoprojectionconverter.o -llas -o $@ ${LIBS} ${LASLIBS} $(INCLUDE) $(LASINCLUDE)
	mkdir -p ../bin	
	cp $@ ../bin
	cp $@ p_las2las
	cp $@ ../bin/p_las2las

laszip: laszip.o geoprojectionconverter.o 
	${LINKER} ${BITS} ${COPTS} laszip.o geoprojectionconverter.o -llas -o $@ ${LIBS} ${LASLIBS} $(INCLUDE) $(LASINCLUDE)
	cp $@ ../bin

lasinfo: lasinfo.o geoprojectionconverter.o
	${LINKER} ${BITS} ${COPTS} lasinfo.o geoprojectionconverter.o -llas -o $@ ${LIBS} ${LASLIBS} $(INCLUDE) $(LASINCLUDE)
	mkdir -p ../bin
	cp $@ ../bin

lasprecision: lasprecision.o geoprojectionconverter.o
	${LINKER} ${BITS} ${COPTS} lasprecision.o geoprojectionconverter.o -llas -o $@ ${LIBS} ${LASLIBS} $(INCLUDE) $(LASINCLUDE)
	mkdir -p ../bin
	cp $@ ../bin

las2txt: las2txt.o
	${LINKER} ${BITS} ${COPTS} las2txt.o -llas -o $@ ${LIBS} ${LASLIBS} $(INCLUDE) $(LASINCLUDE)
	mkdir -p ../bin
	cp $@ ../bin

txt2las: txt2las.o geoprojectionconverter.o
	${LINKER} ${BITS} ${COPTS} txt2las.o geoprojectionconverter.o -llas -o $@ ${LIBS} ${LASLIBS} $(INCLUDE) $(LASINCLUDE)
	mkdir -p ../bin
	cp $@ ../bin

las2las: las2las.o geoprojectionconverter.o
	${LINKER} ${BITS} ${COPTS} las2las.o geoprojectionconverter.o -llas -o $@ ${LIBS} ${LASLIBS} $(INCLUDE) $(LASINCLUDE)
	mkdir -p ../bin
	cp $@ ../bin

lasdiff: lasdiff.o
	${LINKER} ${BITS} ${COPTS} lasdiff.o -llas -o $@ ${LIBS} ${LASLIBS} $(INCLUDE) $(LASINCLUDE)
	mkdir -p ../bin
	cp $@ ../bin

lasmerge: lasmerge.o geoprojectionconverter.o 
	${LINKER} ${BITS} ${COPTS} lasmerge.o geoprojectionconverter.o -llas -o $@ ${LIBS} ${LASLIBS} $(INCLUDE) $(LASINCLUDE)
	mkdir -p ../bin
	cp $@ ../bin

lasindex: lasindex.o
	${LINKER} ${BITS} ${COPTS} lasindex.o -llas -o $@ ${LIBS} ${LASLIBS} $(INCLUDE) $(LASINCLUDE)
	mkdir -p ../bin
	cp $@ ../bin

.cpp.o: 
	${COMPILER} ${BITS} -c ${COPTS} ${INCLUDE} $(LASINCLUDE) $< -o $@

.c.o: 
	${COMPILER} ${BITS} -c ${COPTS} ${INCLUDE} $(LASINCLUDE) $< -o $@

clean:
	rm -rf *.o
	rm -rf laszip lasinfo lasprecision las2txt txt2las las2las lasdiff lasmerge lasindex plas p_las2las

clobber:
	rm -rf *.o
	rm -rf laszip lasinfo lasprecision las2txt txt2las las2las lasdiff lasmerge lasindex plas p_las2las

