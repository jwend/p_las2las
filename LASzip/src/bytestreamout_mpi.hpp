 /*
===============================================================================

  FILE:  bytestreamout_mpi.hpp
  
  CONTENTS:
      
    Class for MPI_File based output streams with endian handling.

  PROGRAMMERS:

    martin.isenburg@rapidlasso.com  -  http://rapidlasso.com
    jwendel@usgs.gov - Jeff Wendel, USGS

  COPYRIGHT:

    (c) 2007-2013, martin isenburg, rapidlasso - fast tools to catch reality

    This is free software; you can redistribute and/or modify it under the
    terms of the GNU Lesser General Licence as published by the Free Software
    Foundation. See the COPYING file for more information.

    This software is distributed WITHOUT ANY WARRANTY and without even the
    implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
  
  CHANGE HISTORY:
  
     1 October 2011 -- added 64 bit file support in MSVC 6.0 at McCafe at Hbf Linz
    10 January 2011 -- licensing change for LGPL release and liblas integration
    12 December 2010 -- created from ByteStreamOutMPI after Howard got pushy (-;
  

===============================================================================
*/
#ifndef BYTE_STREAM_OUT_MPI_H
#define BYTE_STREAM_OUT_MPI_H

#include "bytestreamout.hpp"

#include <stdio.h>

#include <mpi.h>

#if defined(_MSC_VER) && (_MSC_VER < 1300)
extern "C" int _cdecl _fseeki64(FILE*, __int64, int);
extern "C" __int64 _cdecl _ftelli64(FILE*);
#endif

class ByteStreamOutMPI : public ByteStreamOut
{
public:
  ByteStreamOutMPI(MPI_File fh);
/* replace a closed FILE* with a reopened FILE* in "ab" mode */
  BOOL refile(MPI_File fh);
/* write a single byte                                       */
  BOOL putByte(U8 byte);
/* write an array of bytes                                   */
  BOOL putBytes(const U8* bytes, U32 num_bytes);
/* is the stream seekable (e.g. standard out is not)         */
  BOOL isSeekable() const;
/* get current position of stream                            */
  I64 tell() const;
/* seek to this position in the stream                       */
  BOOL seek(const I64 position);
/* seek to the end of the file                               */
  BOOL seekEnd();
/* destructor                                                */
  ~ByteStreamOutMPI(){};
protected:
  MPI_File fh;
};

class ByteStreamOutMPILE : public ByteStreamOutMPI
{
public:
  ByteStreamOutMPILE(MPI_File fh);
/* write 16 bit low-endian field                             */
  BOOL put16bitsLE(const U8* bytes);
/* write 32 bit low-endian field                             */
  BOOL put32bitsLE(const U8* bytes);
/* write 64 bit low-endian field                             */
  BOOL put64bitsLE(const U8* bytes);
/* write 16 bit big-endian field                             */
  BOOL put16bitsBE(const U8* bytes);
/* write 32 bit big-endian field                             */
  BOOL put32bitsBE(const U8* bytes);
/* write 64 bit big-endian field                             */
  BOOL put64bitsBE(const U8* bytes);
private:
  U8 swapped[8];
};

class ByteStreamOutMPIBE : public ByteStreamOutMPI
{
public:
  ByteStreamOutMPIBE(MPI_File fh);
/* write 16 bit low-endian field                             */
  BOOL put16bitsLE(const U8* bytes);
/* write 32 bit low-endian field                             */
  BOOL put32bitsLE(const U8* bytes);
/* write 64 bit low-endian field                             */
  BOOL put64bitsLE(const U8* bytes);
/* write 16 bit big-endian field                             */
  BOOL put16bitsBE(const U8* bytes);
/* write 32 bit big-endian field                             */
  BOOL put32bitsBE(const U8* bytes);
/* write 64 bit big-endian field                             */
  BOOL put64bitsBE(const U8* bytes);
private:
  U8 swapped[8];
};

inline ByteStreamOutMPI::ByteStreamOutMPI(MPI_File fh)
{
  this->fh = fh;
}

inline BOOL ByteStreamOutMPI::refile(MPI_File fh)
{
  if (fh == 0) return FALSE;
  this->fh = fh;
  return TRUE;
}

inline BOOL ByteStreamOutMPI::putByte(U8 byte)
{
  //return (fputc(byte, file) == byte);
  MPI_Status status;
  int rtn = MPI_File_write(fh, &byte, 1, MPI_BYTE, &status);
  if (rtn == MPI_SUCCESS)
    return TRUE;
  else
    return FALSE;

}

inline BOOL ByteStreamOutMPI::putBytes(const U8* bytes, U32 num_bytes)
{
  //return (fwrite(bytes, 1, num_bytes, file) == num_bytes);
  MPI_Status status;
  int rtn = MPI_File_write(fh, (void *)bytes, num_bytes, MPI_BYTE, &status);
  if (rtn == MPI_SUCCESS)
    return TRUE;
  else
    return FALSE;
}

inline BOOL ByteStreamOutMPI::isSeekable() const
{
  //return (file != stdout);
  return TRUE;
}

inline I64 ByteStreamOutMPI::tell() const
{
//#if defined _WIN32 && ! defined (__MINGW32__)
//  return _ftelli64(file);
//#elif defined (__MINGW32__)
//  return (I64)ftello64(file);
//#else
//  return (I64)ftello(file);
//#endif
  MPI_Offset offset;
  MPI_File_get_position(fh, &offset);
  return (I64)offset;


}

inline BOOL ByteStreamOutMPI::seek(I64 position)
{
//#if defined _WIN32 && ! defined (__MINGW32__)
//  return !(_fseeki64(file, position, SEEK_SET));
//#elif defined (__MINGW32__)
//  return !(fseeko64(file, (off_t)position, SEEK_SET));
//#else
//  return !(fseeko(file, (off_t)position, SEEK_SET));
//#endif

  int rtn = MPI_File_seek(fh, position, MPI_SEEK_SET);
  if (rtn == MPI_SUCCESS)
    return TRUE;
  else
    return FALSE;
}

inline BOOL ByteStreamOutMPI::seekEnd()
{
//#if defined _WIN32 && ! defined (__MINGW32__)
//  return !(_fseeki64(file, 0, SEEK_END));
//#elif defined (__MINGW32__)
//  return !(fseeko64(file, (off_t)0, SEEK_END));
//#else
//  return !(fseeko(file, (off_t)0, SEEK_END));
//#endif

  int rtn = MPI_File_seek(fh, 0, MPI_SEEK_END);
  if (rtn == MPI_SUCCESS)
    return TRUE;
  else
    return FALSE;


}

inline ByteStreamOutMPILE::ByteStreamOutMPILE(MPI_File fh) : ByteStreamOutMPI(fh)
{
}

inline BOOL ByteStreamOutMPILE::put16bitsLE(const U8* bytes)
{
  return putBytes(bytes, 2);
}

inline BOOL ByteStreamOutMPILE::put32bitsLE(const U8* bytes)
{
  return putBytes(bytes, 4);
}

inline BOOL ByteStreamOutMPILE::put64bitsLE(const U8* bytes)
{
  return putBytes(bytes, 8);
}

inline BOOL ByteStreamOutMPILE::put16bitsBE(const U8* bytes)
{
  swapped[0] = bytes[1];
  swapped[1] = bytes[0];
  return putBytes(swapped, 2);
}

inline BOOL ByteStreamOutMPILE::put32bitsBE(const U8* bytes)
{
  swapped[0] = bytes[3];
  swapped[1] = bytes[2];
  swapped[2] = bytes[1];
  swapped[3] = bytes[0];
  return putBytes(swapped, 4);
}

inline BOOL ByteStreamOutMPILE::put64bitsBE(const U8* bytes)
{
  swapped[0] = bytes[7];
  swapped[1] = bytes[6];
  swapped[2] = bytes[5];
  swapped[3] = bytes[4];
  swapped[4] = bytes[3];
  swapped[5] = bytes[2];
  swapped[6] = bytes[1];
  swapped[7] = bytes[0];
  return putBytes(swapped, 8);
}

inline ByteStreamOutMPIBE::ByteStreamOutMPIBE(MPI_File fh) : ByteStreamOutMPI(fh)
{
}

inline BOOL ByteStreamOutMPIBE::put16bitsLE(const U8* bytes)
{
  swapped[0] = bytes[1];
  swapped[1] = bytes[0];
  return putBytes(swapped, 2);
}

inline BOOL ByteStreamOutMPIBE::put32bitsLE(const U8* bytes)
{
  swapped[0] = bytes[3];
  swapped[1] = bytes[2];
  swapped[2] = bytes[1];
  swapped[3] = bytes[0];
  return putBytes(swapped, 4);
}

inline BOOL ByteStreamOutMPIBE::put64bitsLE(const U8* bytes)
{
  swapped[0] = bytes[7];
  swapped[1] = bytes[6];
  swapped[2] = bytes[5];
  swapped[3] = bytes[4];
  swapped[4] = bytes[3];
  swapped[5] = bytes[2];
  swapped[6] = bytes[1];
  swapped[7] = bytes[0];
  return putBytes(swapped, 8);
}

inline BOOL ByteStreamOutMPIBE::put16bitsBE(const U8* bytes)
{
  return putBytes(bytes, 2);
}

inline BOOL ByteStreamOutMPIBE::put32bitsBE(const U8* bytes)
{
  return putBytes(bytes, 4);
}

inline BOOL ByteStreamOutMPIBE::put64bitsBE(const U8* bytes)
{
  return putBytes(bytes, 8);
}

#endif
