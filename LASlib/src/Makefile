# makefile for liblas.a
#
#COPTS    = -g -Wall
COPTS     = -O3 -Wall -Wno-deprecated -Wno-strict-aliasing -Wno-unused-but-set-variable -Wno-unused-result -DNDEBUG -DUNORDERED
#COPTS  = -O3 -Wno-deprecated
COMPILER  ?= mpic++
AR  = ar rc
#BITS     = -64

INCLUDE		= -I/usr/include/ -I../zip/src -I../inc -I.

OBJ_LAS		= lasreader.o laswriter.o lasreader_las.o lasreader_bin.o lasreader_qfit.o lasreader_shp.o lasreader_asc.o lasreader_bil.o lasreader_dtm.o lasreader_txt.o lasreadermerged.o lasreaderbuffered.o lasreaderpipeon.o laswriter_las.o laswriter_bin.o laswriter_qfit.o laswriter_wrl.o laswriter_txt.o laswaveform13reader.o laswaveform13writer.o lasutility.o lasinterval.o lasspatial.o lasquadtree.o lasindex.o lasfilter.o lastransform.o fopen_compressed.o

OBJ_LAZ		= ../zip/src/laszip.o ../zip/src/lasreadpoint.o ../zip/src/lasreaditemcompressed_v1.o ../zip/src/lasreaditemcompressed_v2.o ../zip/src/laswritepoint.o  ../zip/src/laswriteitemcompressed_v1.o ../zip/src/laswriteitemcompressed_v2.o ../zip/src/integercompressor.o ../zip/src/arithmeticdecoder.o ../zip/src/arithmeticencoder.o ../zip/src/arithmeticmodel.o

all: liblas.a

liblas.a: ${OBJ_LAZ} ${OBJ_LAS} 
	$(AR) $@ ${OBJ_LAZ} ${OBJ_LAS}
	cp -p $@ ../lib

.cpp.o:
	${COMPILER} ${BITS} -c ${COPTS} ${INCLUDE} $< -o $@

.c.o: 
	${COMPILER} ${BITS} -c ${COPTS} ${INCLUDE} $< -o $@

clean:
	rm -rf ${OBJ_LAZ} ${OBJ_LAS} 

clobber:
	rm -rf liblas.a
	rm -rf ${OBJ_LAZ} ${OBJ_LAS}

