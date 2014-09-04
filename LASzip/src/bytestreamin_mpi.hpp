/*
===============================================================================

  FILE:  bytestreamin_mpi.hpp
  
  CONTENTS:
      
    Class for MPI_Fileinput streams with endian handling.

  PROGRAMMERS:

    martin.isenburg@rapidlasso.com  -  http://rapidlasso.com
    jwendel@usgs.gov - Jeff Wendel, USGS
  COPYRIGHT:

    (c) 2007-2012, martin isenburg, rapidlasso - tools to catch reality

    This is free software; you can redistribute and/or modify it under the
    terms of the GNU Lesser General Licence as published by the Free Software
    Foundation. See the COPYING file for more information.

    This software is distributed WITHOUT ANY WARRANTY and without even the
    implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
  
  CHANGE HISTORY:
  
     1 October 2011 -- added 64 bit file support in MSVC 6.0 at McCafe at Hbf Linz
    10 January 2011 -- licensing change for LGPL release and liblas integration
    12 December 2010 -- created from ByteStreamOutFile after Howard got pushy (-;
  
===============================================================================
*/
#ifndef BYTE_STREAM_IN_MPI_H
#define BYTE_STREAM_IN_MPI_H

#include "bytestreamin.hpp"

#include <stdio.h>
#include <mpi.h>

#if defined(_MSC_VER) && (_MSC_VER < 1300)
extern "C" __int64 _cdecl _ftelli64(FILE*);
extern "C" int _cdecl _fseeki64(FILE*, __int64, int);
#endif

class ByteStreamInMPI : public ByteStreamIn
{
public:
  ByteStreamInMPI(MPI_File fh);
/* read a single byte                                        */
  U32 getByte();
/* read an array of bytes                                    */
  void getBytes(U8* bytes, const U32 num_bytes);
/* is the stream seekable (e.g. stdin is not)                */
  BOOL isSeekable() const;
/* get current position of stream                            */
  I64 tell() const;
/* seek to this position in the stream                       */
  BOOL seek(const I64 position);
/* seek to the end of the file                               */
  BOOL seekEnd(const I64 distance=0);
/* destructor                                                */
  ~ByteStreamInMPI(){};
protected:
  MPI_File fh;
};

class ByteStreamInMPILE : public ByteStreamInMPI
{
public:
  ByteStreamInMPILE(MPI_File fh);
/* read 16 bit low-endian field                              */
  void get16bitsLE(U8* bytes);
/* read 32 bit low-endian field                              */
  void get32bitsLE(U8* bytes);
/* read 64 bit low-endian field                              */
  void get64bitsLE(U8* bytes);
/* read 16 bit big-endian field                              */
  void get16bitsBE(U8* bytes);
/* read 32 bit big-endian field                              */
  void get32bitsBE(U8* bytes);
/* read 64 bit big-endian field                              */
  void get64bitsBE(U8* bytes);
private:
  U8 swapped[8];
};

class ByteStreamInMPIBE : public ByteStreamInMPI
{
public:
  ByteStreamInMPIBE(MPI_File fh);
/* read 16 bit low-endian field                              */
  void get16bitsLE(U8* bytes);
/* read 32 bit low-endian field                              */
  void get32bitsLE(U8* bytes);
/* read 64 bit low-endian field                              */
  void get64bitsLE(U8* bytes);
/* read 16 bit big-endian field                              */
  void get16bitsBE(U8* bytes);
/* read 32 bit big-endian field                              */
  void get32bitsBE(U8* bytes);
/* read 64 bit big-endian field                              */
  void get64bitsBE(U8* bytes);
private:
  U8 swapped[8];
};

inline ByteStreamInMPI::ByteStreamInMPI(MPI_File fh)
{
  this->fh = fh;
}

inline U32 ByteStreamInMPI::getByte()
{
  //int byte = getc(file);
  //if (byte == EOF)
  //{
  //  throw EOF;
  //}
  //return (U32)byte;
  U32 byte;
  MPI_Status status;
  int cnt;
  MPI_File_read(fh, &byte, 1, MPI_BYTE, &status);
  MPI_Get_count(&status, MPI_BYTE, &cnt);
  if(cnt != 1){
      throw EOF;
  }
  return byte;
}

inline void ByteStreamInMPI::getBytes(U8* bytes, const U32 num_bytes)
{
//  if (fread(bytes, 1, num_bytes, file) != num_bytes)
//  {
//    throw EOF;
//  }

	MPI_Status status;
	int cnt;
	MPI_File_read(fh, bytes, num_bytes, MPI_BYTE, &status);
	MPI_Get_count(&status, MPI_BYTE, &cnt);
	if(cnt != num_bytes){
		throw EOF;
	}
}

inline BOOL ByteStreamInMPI::isSeekable() const
{
  return TRUE;
}

inline I64 ByteStreamInMPI::tell() const
{
	  MPI_Offset offset;
	  MPI_File_get_position(fh, &offset);
	  return (I64)offset;
}

inline BOOL ByteStreamInMPI::seek(const I64 position)
{

  int rtn = MPI_File_seek(fh, position, MPI_SEEK_SET);
  if (rtn == MPI_SUCCESS)
    return TRUE;
  else
    return FALSE;

}

inline BOOL ByteStreamInMPI::seekEnd(const I64 distance)
{
	 int rtn = MPI_File_seek(fh, 0, MPI_SEEK_END);
	  if (rtn == MPI_SUCCESS)
	    return TRUE;
	  else
	    return FALSE;

}

inline ByteStreamInMPILE::ByteStreamInMPILE(MPI_File fh) : ByteStreamInMPI(fh)
{
}

inline void ByteStreamInMPILE::get16bitsLE(U8* bytes)
{
  getBytes(bytes, 2);
}

inline void ByteStreamInMPILE::get32bitsLE(U8* bytes)
{
  getBytes(bytes, 4);
}

inline void ByteStreamInMPILE::get64bitsLE(U8* bytes)
{
  getBytes(bytes, 8);
}

inline void ByteStreamInMPILE::get16bitsBE(U8* bytes)
{
  getBytes(swapped, 2);
  bytes[0] = swapped[1];
  bytes[1] = swapped[0];
}

inline void ByteStreamInMPILE::get32bitsBE(U8* bytes)
{
  getBytes(swapped, 4);
  bytes[0] = swapped[3];
  bytes[1] = swapped[2];
  bytes[2] = swapped[1];
  bytes[3] = swapped[0];
}

inline void ByteStreamInMPILE::get64bitsBE(U8* bytes)
{
  getBytes(swapped, 8);
  bytes[0] = swapped[7];
  bytes[1] = swapped[6];
  bytes[2] = swapped[5];
  bytes[3] = swapped[4];
  bytes[4] = swapped[3];
  bytes[5] = swapped[2];
  bytes[6] = swapped[1];
  bytes[7] = swapped[0];
}

inline ByteStreamInMPIBE::ByteStreamInMPIBE(MPI_File fh) : ByteStreamInMPI(fh)
{
}

inline void ByteStreamInMPIBE::get16bitsLE(U8* bytes)
{
  getBytes(swapped, 2);
  bytes[0] = swapped[1];
  bytes[1] = swapped[0];
}

inline void ByteStreamInMPIBE::get32bitsLE(U8* bytes)
{
  getBytes(swapped, 4);
  bytes[0] = swapped[3];
  bytes[1] = swapped[2];
  bytes[2] = swapped[1];
  bytes[3] = swapped[0];
}

inline void ByteStreamInMPIBE::get64bitsLE(U8* bytes)
{
  getBytes(swapped, 8);
  bytes[0] = swapped[7];
  bytes[1] = swapped[6];
  bytes[2] = swapped[5];
  bytes[3] = swapped[4];
  bytes[4] = swapped[3];
  bytes[5] = swapped[2];
  bytes[6] = swapped[1];
  bytes[7] = swapped[0];
}

inline void ByteStreamInMPIBE::get16bitsBE(U8* bytes)
{
  getBytes(bytes, 2);
}

inline void ByteStreamInMPIBE::get32bitsBE(U8* bytes)
{
  getBytes(bytes, 4);
}

inline void ByteStreamInMPIBE::get64bitsBE(U8* bytes)
{
  getBytes(bytes, 8);
}

#endif
