/*
===============================================================================

  FILE:  las2las.cpp
  
  CONTENTS:
  
    This tool reads and writes LIDAR data in the LAS format and is typically
    used to modify the contents of a LAS file. Examples are keeping or dropping
    those points lying within a certain region (-keep_xy, -drop_x_above, ...),
    or points with a certain elevation (-keep_z, -drop_z_below, -drop_z_above)
    eliminating points that are the second return (-drop_return 2), that have a
    scan angle above a certain threshold (-drop_scan_angle_above 5), or that are
    below a certain intensity (-drop_intensity_below 15).
    Another typical use may be to extract only first (-first_only) returns or
    only last returns (-last_only). Extracting the first return is actually the
    same as eliminating all others (e.g. -keep_return 2 -keep_return 3, ...).

  PROGRAMMERS:
  
    martin.isenburg@rapidlasso.com  -  http://rapidlasso.com
  
  COPYRIGHT:
  
    (c) 2007-14, martin isenburg, rapidlasso - fast tools to catch reality

    This is free software; you can redistribute and/or modify it under the
    terms of the GNU Lesser General Licence as published by the Free Software
    Foundation. See the LICENSE.txt file for more information.

    This software is distributed WITHOUT ANY WARRANTY and without even the
    implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
  
  CHANGE HISTORY:
  
     5 July 2012 -- added option to '-remove_original_vlr' 
     6 May 2012 -- added option to '-remove_tiling_vlr' 
     5 January 2012 -- added option to clip points to the bounding box
    17 May 2011 -- enabling batch processing with wildcards or multiple file names
    13 May 2011 -- moved indexing, filtering, transforming into LASreader
    18 April 2011 -- can set projection tags or reproject horizontally
    26 January 2011 -- added LAStransform for simply manipulations of points 
    21 January 2011 -- added LASreadOpener and reading of multiple LAS files 
     3 January 2011 -- added -reoffset & -rescale + -keep_circle via LASfilter
    10 January 2010 -- added -subseq for selecting a [start, end] interval
    10 June 2009 -- added -scale_rgb_down and -scale_rgb_up to scale rgb values
    12 March 2009 -- updated to ask for input if started without arguments 
     9 March 2009 -- added ability to remove user defined headers or vlrs
    17 September 2008 -- updated to deal with LAS format version 1.2
    17 September 2008 -- dropping or keeping in double precision and based on z
    10 July 2007 -- created after talking with Linda about the H1B process
  
===============================================================================
*/

#include <time.h>
#include <stdlib.h>
#include <string.h>
#include <mpi.h>
#include <math.h>
#include "lasreader.hpp"
#include "bytestreamin.hpp"
#include "laswriter.hpp"
#include "geoprojectionconverter.hpp"

static void usage(bool error=false, bool wait=false)
{
  fprintf(stderr,"usage:\n");
  fprintf(stderr,"las2las -i *.las -utm 13N\n");
  fprintf(stderr,"las2las -i *.laz -first_only -olaz\n");
  fprintf(stderr,"las2las -i *.las -drop_return 4 5 -olaz\n");
  fprintf(stderr,"las2las -latlong -target_utm 12T -i in.las -o out.las\n");
  fprintf(stderr,"las2las -point_type 0 -lof file_list.txt -merged -o out.las\n");
  fprintf(stderr,"las2las -remove_vlr 2 -scale_rgb_up -i in.las -o out.las\n");
  fprintf(stderr,"las2las -i in.las -keep_xy 630000 4834500 630500 4835000 -keep_z 10 100 -o out.las\n");
  fprintf(stderr,"las2las -i in.txt -iparse xyzit -keep_circle 630200 4834750 100 -oparse xyzit -o out.txt\n");
  fprintf(stderr,"las2las -i in.las -keep_scan_angle -15 15 -o out.las\n");
  fprintf(stderr,"las2las -i in.las -rescale 0.01 0.01 0.01 -reoffset 0 300000 0 -o out.las\n");
  fprintf(stderr,"las2las -i in.las -set_version 1.2 -keep_gpstime 46.5 47.5 -o out.las\n");
  fprintf(stderr,"las2las -i in.las -drop_intensity_below 10 -olaz -stdout > out.laz\n");
  fprintf(stderr,"las2las -i in.las -last_only -drop_gpstime_below 46.75 -otxt -oparse xyzt -stdout > out.txt\n");
  fprintf(stderr,"las2las -i in.las -remove_all_vlrs -keep_class 2 3 4 -olas -stdout > out.las\n");
  fprintf(stderr,"las2las -h\n");
  if (wait)
  {
    fprintf(stderr,"<press ENTER>\n");
    getc(stdin);
  }
  exit(error);
}

static void byebye(bool error=false, bool wait=false)
{
  if (wait)
  {
    fprintf(stderr,"<press ENTER>\n");
    getc(stdin);
  }
  exit(error);
}

static double taketime()
{
  return (double)(clock())/CLOCKS_PER_SEC;
}

// for point type conversions
const U8 convert_point_type_from_to[11][11] = 
{
  {  0,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1 },
  {  0,  0,  1,  1,  1,  1,  1,  1,  1,  1,  1 },
  {  0,  1,  0,  1,  1,  1,  1,  1,  1,  1,  1 },
  {  0,  0,  1,  0,  1,  1,  1,  1,  1,  1,  1 },
  {  0,  0,  1,  1,  0,  1,  1,  1,  1,  1,  1 },
  {  0,  0,  1,  0,  1,  0,  1,  1,  1,  1,  1 },
  {  1,  1,  1,  1,  1,  1,  0,  1,  1,  1,  1 },
  {  1,  1,  1,  1,  1,  1,  1,  0,  1,  1,  1 },
  {  1,  1,  1,  1,  1,  1,  1,  1,  0,  1,  1 },
  {  1,  1,  1,  1,  1,  1,  1,  1,  1,  0,  1 },
  {  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  0 },
};


int main(int argc, char *argv[])
{
  int i;
  int is_mpi = 1;
  int debug = 1;
  bool verbose = false;
  bool force = false;
  // fixed header changes 
  int set_version_major = -1;
  int set_version_minor = -1;
  int set_point_data_format = -1;
  int set_point_data_record_length = -1;
  int set_gps_time_endcoding = -1;
  // variable header changes
  bool remove_extra_header = false;
  bool remove_all_variable_length_records = false;
  int remove_variable_length_record = -1;
  int remove_variable_length_record_from = -1;
  int remove_variable_length_record_to = -1;
  bool remove_tiling_vlr = false;
  bool remove_original_vlr = false;
  // extract a subsequence
  //unsigned int subsequence_start = 0;
  //unsigned int subsequence_stop = U32_MAX;
  I64 subsequence_start = 0;
  I64 subsequence_stop = I64_MAX;


  // fix files with corrupt points
  bool clip_to_bounding_box = false;
  double start_time = 0;
  time_t wall_start_time;
  time_t wall_end_time;
  LASreadOpener lasreadopener;
  //if(is_mpi)lasreadopener.setIsMpi(TRUE);
  GeoProjectionConverter geoprojectionconverter;
  LASwriteOpener laswriteopener;
  if(is_mpi)laswriteopener.setIsMpi(TRUE);


  int process_count = 1;
  int rank = 0;
  start_time = taketime();
  time(&wall_start_time);

  if (is_mpi){
      MPI_Init(&argc,&argv);
      MPI_Comm_size(MPI_COMM_WORLD,&process_count);
      MPI_Comm_rank(MPI_COMM_WORLD,&rank);
      if(debug) printf ("MPI task %d has started...\n", rank);
  }



  if (argc == 1)
  {

    fprintf(stderr,"las2las.exe is better run in the command line or via the lastool.exe GUI\n");
    char file_name[256];
    fprintf(stderr,"enter input file: "); fgets(file_name, 256, stdin);
    file_name[strlen(file_name)-1] = '\0';
    lasreadopener.set_file_name(file_name);
    fprintf(stderr,"enter output file: "); fgets(file_name, 256, stdin);
    file_name[strlen(file_name)-1] = '\0';
    laswriteopener.set_file_name(file_name);

  }
  else
  {
    for (i = 1; i < argc; i++)
    {
      //if (argv[i][0] == '�') argv[i][0] = '-';
      if (strcmp(argv[i],"-week_to_adjusted") == 0)
      {
        set_gps_time_endcoding = 1;
      }
      else if (strcmp(argv[i],"-adjusted_to_week") == 0)
      {
        set_gps_time_endcoding = 0;
      }
    }
    if (!geoprojectionconverter.parse(argc, argv)) byebye(true);
    if (!lasreadopener.parse(argc, argv)) byebye(true);
    if (!laswriteopener.parse(argc, argv)) byebye(true);
  }

  for (i = 1; i < argc; i++)
  {
    if (argv[i][0] == '\0')
    {
      continue;
    }
    else if (strcmp(argv[i],"-h") == 0 || strcmp(argv[i],"-help") == 0)
    {
      fprintf(stderr, "LAStools (by martin@rapidlasso.com) version %d\n", LAS_TOOLS_VERSION);
      usage();
    }
    else if (strcmp(argv[i],"-v") == 0 || strcmp(argv[i],"-verbose") == 0)
    {
      verbose = true;
    }
    else if (strcmp(argv[i],"-version") == 0)
    {
      fprintf(stderr, "LAStools (by martin@rapidlasso.com) version %d\n", LAS_TOOLS_VERSION);
      byebye();
    }
    else if (strcmp(argv[i],"-gui") == 0)
    {

      fprintf(stderr, "WARNING: not compiled with GUI support. ignoring '-gui' ...\n");

    }
    else if (strcmp(argv[i],"-cores") == 0)
    {

      fprintf(stderr, "WARNING: not compiled with multi-core batching. ignoring '-cores' ...\n");
      i++;

    }
    else if (strcmp(argv[i],"-force") == 0)
    {
      force = true;
    }
    else if (strcmp(argv[i],"-subseq") == 0)
    {
      if ((i+2) >= argc)
      {
        fprintf(stderr,"ERROR: '%s' needs 2 arguments: start stop\n", argv[i]);
        byebye(true);
      }
      subsequence_start = (unsigned int)atoi(argv[i+1]); subsequence_stop = (unsigned int)atoi(argv[i+2]);
      i+=2;
    }
    else if (strcmp(argv[i],"-start_at_point") == 0)
    {
      if ((i+1) >= argc)
      {
        fprintf(stderr,"ERROR: '%s' needs 1 argument: start\n", argv[i]);
        byebye(true);
      }
      subsequence_start = (unsigned int)atoi(argv[i+1]);
      i+=1;
    }
    else if (strcmp(argv[i],"-stop_at_point") == 0)
    {
      if ((i+1) >= argc)
      {
        fprintf(stderr,"ERROR: '%s' needs 1 argument: stop\n", argv[i]);
        byebye(true);
      }
      subsequence_stop = (unsigned int)atoi(argv[i+1]);
      i+=1;
    }
    else if (strcmp(argv[i],"-set_version") == 0)
    {
      if ((i+1) >= argc)
      {
        fprintf(stderr,"ERROR: '%s' needs 1 argument: major.minor\n", argv[i]);
        byebye(true);
      }
      if (sscanf(argv[i+1],"%d.%d",&set_version_major,&set_version_minor) != 2)
      {
        fprintf(stderr, "ERROR: cannot understand argument '%s' for '%s'\n", argv[i+1], argv[i]);
        usage(true);
      }
      i+=1;
    }
    else if (strcmp(argv[i],"-set_version_major") == 0)
    {
      if ((i+1) >= argc)
      {
        fprintf(stderr,"ERROR: '%s' needs 1 argument: major\n", argv[i]);
        byebye(true);
      }
      set_version_major = atoi(argv[i+1]);
      i+=1;
    }
    else if (strcmp(argv[i],"-set_version_minor") == 0)
    {
      if ((i+1) >= argc)
      {
        fprintf(stderr,"ERROR: '%s' needs 1 argument: minor\n", argv[i]);
        byebye(true);
      }
      set_version_minor = atoi(argv[i+1]);
      i+=1;
    }
    else if (strcmp(argv[i],"-remove_extra") == 0)
    {
      remove_extra_header = true;
    }
    else if (strcmp(argv[i],"-remove_all_vlrs") == 0)
    {
      remove_all_variable_length_records = true;
    }
    else if (strcmp(argv[i],"-remove_vlr") == 0)
    {
      if ((i+1) >= argc)
      {
        fprintf(stderr,"ERROR: '%s' needs 1 argument: number\n", argv[i]);
        byebye(true);
      }
      remove_variable_length_record = atoi(argv[i+1]);
      remove_variable_length_record_from = -1;
      remove_variable_length_record_to = -1;
      i++;
    }
    else if (strcmp(argv[i],"-remove_vlrs_from_to") == 0)
    {
      if ((i+2) >= argc)
      {
        fprintf(stderr,"ERROR: '%s' needs 2 arguments: start end\n", argv[i]);
        byebye(true);
      }
      remove_variable_length_record = -1;
      remove_variable_length_record_from = atoi(argv[i+1]);
      remove_variable_length_record_to = atoi(argv[i+2]);
      i+=2;
    }
    else if (strcmp(argv[i],"-remove_tiling_vlr") == 0)
    {
      remove_tiling_vlr = true;
      i++;
    }
    else if (strcmp(argv[i],"-remove_original_vlr") == 0)
    {
      remove_original_vlr = true;
      i++;
    }
    else if (strcmp(argv[i],"-set_point_type") == 0 || strcmp(argv[i],"-set_point_data_format") == 0 || strcmp(argv[i],"-point_type") == 0) 
    {
      if ((i+1) >= argc)
      {
        fprintf(stderr,"ERROR: '%s' needs 1 argument: type\n", argv[i]);
        byebye(true);
      }
      set_point_data_format = atoi(argv[i+1]);
      i++;
    }
    else if (strcmp(argv[i],"-set_point_data_record_length") == 0 || strcmp(argv[i],"-set_point_size") == 0 || strcmp(argv[i],"-point_size") == 0) 
    {
      if ((i+1) >= argc)
      {
        fprintf(stderr,"ERROR: '%s' needs 1 argument: size\n", argv[i]);
        byebye(true);
      }
      set_point_data_record_length = atoi(argv[i+1]);
      i++;
    }
    else if (strcmp(argv[i],"-clip_to_bounding_box") == 0 || strcmp(argv[i],"-clip_to_bb") == 0) 
    {
      clip_to_bounding_box = true;
    }
    else if ((argv[i][0] != '-') && (lasreadopener.get_file_name_number() == 0))
    {
      lasreadopener.add_file_name(argv[i]);
      argv[i][0] = '\0';
    }
    else
    {
      fprintf(stderr, "ERROR: cannot understand argument '%s'\n", argv[i]);
      usage(true);
    }
  }



  // check input

  if (!lasreadopener.active())
  {
    fprintf(stderr,"ERROR: no input specified\n");
    usage(true, argc==1);
  }
  
  BOOL extra_pass = laswriteopener.is_piped();

  // for piped output we need an extra pass

  if (extra_pass)
  {
    if (lasreadopener.is_piped())
    {
      fprintf(stderr, "ERROR: input and output cannot both be piped\n");
      usage(true);
    }
  }

  // make sure we do not corrupt the input file

  if (lasreadopener.get_file_name() && laswriteopener.get_file_name() && (strcmp(lasreadopener.get_file_name(), laswriteopener.get_file_name()) == 0))
  {
    fprintf(stderr, "ERROR: input and output file name are identical\n");
    usage(true);
  }
    
  // possibly loop over multiple input files

  while (lasreadopener.active())
  {
   // if (verbose) start_time = taketime();

    // open lasreader

    LASreader* lasreader = lasreadopener.open();

    if (lasreader == 0)
    {
      fprintf(stderr, "ERROR: could not open lasreader\n");
      usage(true, argc==1);
    }

    // store the inventory for the header

    LASinventory lasinventory;

    // the point we write sometimes needs to be copied

    LASpoint* point = 0;

    // prepare the header for output

    if (set_gps_time_endcoding != -1)
    {
      if (set_gps_time_endcoding == 0)
      {
        if ((lasreader->header.global_encoding & 1) == 0)
        {
          fprintf(stderr, "WARNING: global encoding indicates file already in GPS week time\n");
          if (force)
          {
            fprintf(stderr, "         forced conversion.\n");
          }
          else
          {
            fprintf(stderr, "         use '-force' to force conversion.\n");
            byebye(true);
          }
        }
        else
        {
          lasreader->header.global_encoding &= ~1;
        }
      }
      else if (set_gps_time_endcoding == 1)
      {
        if ((lasreader->header.global_encoding & 1) == 1)
        {
          fprintf(stderr, "WARNING: global encoding indicates file already in Adjusted Standard GPS time\n");
          if (force)
          {
            fprintf(stderr, "         forced conversion.\n");
          }
          else
          {
            fprintf(stderr, "         use '-force' to force conversion.\n");
            byebye(true);
          }
        }
        else
        {
          lasreader->header.global_encoding |= 1;
        }
      }
    }

    if (set_version_major != -1)
    {
      if (set_version_major != 1)
      {
        fprintf(stderr, "ERROR: unknown version_major %d\n", set_version_major);
        byebye(true);
      }
      lasreader->header.version_major = (U8)set_version_major;
    }

    if (set_version_minor >= 0)
    {
      if (set_version_minor > 4)
      {
        fprintf(stderr, "ERROR: unknown version_minor %d\n", set_version_minor);
        byebye(true);
      }
      if (set_version_minor < 3)
      {
        if (lasreader->header.version_minor == 3)
        {
          lasreader->header.header_size -= 8;
          lasreader->header.offset_to_point_data -= 8;
        }
        else if (lasreader->header.version_minor >= 4)
        {
          lasreader->header.header_size -= (8 + 140);
          lasreader->header.offset_to_point_data -= (8 + 140);
        }
      }
      else if (set_version_minor == 3)
      {
        if (lasreader->header.version_minor < 3)
        {
          lasreader->header.header_size += 8;
          lasreader->header.offset_to_point_data += 8;
          lasreader->header.start_of_waveform_data_packet_record = 0;
        }
        else if (lasreader->header.version_minor >= 4)
        {
          lasreader->header.header_size -= 140;
          lasreader->header.offset_to_point_data -= 140;
        }
      }
      else if (set_version_minor == 4) 
      {
        if (lasreader->header.version_minor < 3)
        {
          lasreader->header.header_size += (8 + 140);
          lasreader->header.offset_to_point_data += (8 + 140);
          lasreader->header.start_of_waveform_data_packet_record = 0;
        }
        else if (lasreader->header.version_minor == 3)
        {
          lasreader->header.header_size += 140;
          lasreader->header.offset_to_point_data += 140;
        }
      }

      if ((set_version_minor <= 3) && (lasreader->header.version_minor >= 4))
      {
        if (lasreader->header.point_data_format > 5)
        {
          switch (lasreader->header.point_data_format)
          {
          case 6:
            fprintf(stderr, "WARNING: downgrading point_data_format from %d to 1\n", lasreader->header.point_data_format);
            lasreader->header.point_data_format = 1;
            fprintf(stderr, "         and point_data_record_length from %d to %d\n", lasreader->header.point_data_record_length, lasreader->header.point_data_record_length - 2);
            lasreader->header.point_data_record_length -= 2;
            break;
          case 7:
            fprintf(stderr, "WARNING: downgrading point_data_format from %d to 3\n", lasreader->header.point_data_format);
            lasreader->header.point_data_format = 3;
            fprintf(stderr, "         and point_data_record_length from %d to %d\n", lasreader->header.point_data_record_length, lasreader->header.point_data_record_length - 2);
            lasreader->header.point_data_record_length -= 2;
            break;
          case 8:
            fprintf(stderr, "WARNING: downgrading point_data_format from %d to 3\n", lasreader->header.point_data_format);
            lasreader->header.point_data_format = 3;
            fprintf(stderr, "         and point_data_record_length from %d to %d\n", lasreader->header.point_data_record_length, lasreader->header.point_data_record_length - 4);
            lasreader->header.point_data_record_length -= 4;
            break;
          case 9:
            fprintf(stderr, "WARNING: downgrading point_data_format from %d to 4\n", lasreader->header.point_data_format);
            lasreader->header.point_data_format = 4;
            fprintf(stderr, "         and point_data_record_length from %d to %d\n", lasreader->header.point_data_record_length, lasreader->header.point_data_record_length - 2);
            lasreader->header.point_data_record_length -= 2;
            break;
          case 10:
            fprintf(stderr, "WARNING: downgrading point_data_format from %d to 5\n", lasreader->header.point_data_format);
            lasreader->header.point_data_format = 5;
            fprintf(stderr, "         and point_data_record_length from %d to %d\n", lasreader->header.point_data_record_length, lasreader->header.point_data_record_length - 4);
            lasreader->header.point_data_record_length -= 4;
            break;
          default:
            fprintf(stderr, "ERROR: unknown point_data_format %d\n", lasreader->header.point_data_format);
            byebye(true);
          }
        }
        point = new LASpoint;
        point->init(&lasreader->header, lasreader->header.point_data_format, lasreader->header.point_data_record_length);
      }

      lasreader->header.version_minor = (U8)set_version_minor;
    }

    // are we supposed to change the point data format

    if (set_point_data_format != -1)
    {
      if (set_point_data_format < 0 || set_point_data_format > 10)
      {
        fprintf(stderr, "ERROR: unknown point_data_format %d\n", set_point_data_format);
        byebye(true);
      }
      // depending on the conversion we may need to copy the point
      if (convert_point_type_from_to[lasreader->header.point_data_format][set_point_data_format])
      {
        if (point == 0) point = new LASpoint;
      }
      lasreader->header.point_data_format = (U8)set_point_data_format;
      lasreader->header.clean_laszip();
      switch (lasreader->header.point_data_format)
      {
      case 0:
        lasreader->header.point_data_record_length = 20;
        break;
      case 1:
        lasreader->header.point_data_record_length = 28;
        break;
      case 2:
        lasreader->header.point_data_record_length = 26;
        break;
      case 3:
        lasreader->header.point_data_record_length = 34;
        break;
      case 4:
        lasreader->header.point_data_record_length = 57;
        break;
      case 5:
        lasreader->header.point_data_record_length = 63;
        break;
      case 6:
        lasreader->header.point_data_record_length = 30;
        break;
      case 7:
        lasreader->header.point_data_record_length = 36;
        break;
      case 8:
        lasreader->header.point_data_record_length = 38;
        break;
      case 9:
        lasreader->header.point_data_record_length = 59;
        break;
      case 10:
        lasreader->header.point_data_record_length = 67;
        break;
      }
    }

    // are we supposed to change the point data record length

    if (set_point_data_record_length != -1)
    {
      I32 num_extra_bytes = 0;
      switch (lasreader->header.point_data_format)
      {
      case 0:
        num_extra_bytes = set_point_data_record_length - 20;
        break;
      case 1:
        num_extra_bytes = set_point_data_record_length - 28;
        break;
      case 2:
        num_extra_bytes = set_point_data_record_length - 26;
        break;
      case 3:
        num_extra_bytes = set_point_data_record_length - 34;
        break;
      case 4:
        num_extra_bytes = set_point_data_record_length - 57;
        break;
      case 5:
        num_extra_bytes = set_point_data_record_length - 63;
        break;
      case 6:
        num_extra_bytes = set_point_data_record_length - 30;
        break;
      case 7:
        num_extra_bytes = set_point_data_record_length - 36;
        break;
      case 8:
        num_extra_bytes = set_point_data_record_length - 38;
        break;
      case 9:
        num_extra_bytes = set_point_data_record_length - 59;
        break;
      case 10:
        num_extra_bytes = set_point_data_record_length - 67;
        break;
      }
      if (num_extra_bytes < 0)
      {
        fprintf(stderr, "ERROR: point_data_format %d needs record length of at least %d\n", lasreader->header.point_data_format, set_point_data_record_length - num_extra_bytes);
        byebye(true);
      }
      if (lasreader->header.point_data_record_length < set_point_data_record_length)
      {
        if (!point) point = new LASpoint;
      }
      lasreader->header.point_data_record_length = (U16)set_point_data_record_length;
      lasreader->header.clean_laszip();
    }

    // if the point needs to be copied set up the data fields

    if (point)
    {
      point->init(&lasreader->header, lasreader->header.point_data_format, lasreader->header.point_data_record_length);
    }

    // maybe we should remove some stuff

    if (remove_extra_header)
    {
      lasreader->header.clean_user_data_in_header();
      lasreader->header.clean_user_data_after_header();
    }

    if (remove_all_variable_length_records)
    {
      lasreader->header.clean_vlrs();
    }
    else
    {
      if (remove_variable_length_record != -1)
      {
        lasreader->header.remove_vlr(remove_variable_length_record);
      }
    
      if (remove_variable_length_record_from != -1)
      {
        for (i = remove_variable_length_record_to; i >= remove_variable_length_record_from; i--)
        {
          lasreader->header.remove_vlr(i);
        }
      }
    }

    if (remove_tiling_vlr)
    {
      lasreader->header.clean_lastiling();
    }

    if (remove_original_vlr)
    {
      lasreader->header.clean_lasoriginal();
    }

    // maybe we should add / change the projection information
    LASquantizer* reproject_quantizer = 0;
    LASquantizer* saved_quantizer = 0;
    if (geoprojectionconverter.has_projection(true) || geoprojectionconverter.has_projection(false))
    {
      if (!geoprojectionconverter.has_projection(true) && lasreader->header.vlr_geo_keys)
      {
        geoprojectionconverter.set_projection_from_geo_keys(lasreader->header.vlr_geo_keys[0].number_of_keys, (GeoProjectionGeoKeys*)lasreader->header.vlr_geo_key_entries, lasreader->header.vlr_geo_ascii_params, lasreader->header.vlr_geo_double_params);
      }

      if (geoprojectionconverter.has_projection(true) && geoprojectionconverter.has_projection(false))
      {
        reproject_quantizer = new LASquantizer();
        double point[3];
        point[0] = (lasreader->header.min_x+lasreader->header.max_x)/2;
        point[1] = (lasreader->header.min_y+lasreader->header.max_y)/2;
        point[2] = (lasreader->header.min_z+lasreader->header.max_z)/2;
        geoprojectionconverter.to_target(point);
        reproject_quantizer->x_scale_factor = geoprojectionconverter.get_target_precision();
        reproject_quantizer->y_scale_factor = geoprojectionconverter.get_target_precision();
        reproject_quantizer->z_scale_factor = lasreader->header.z_scale_factor;
        reproject_quantizer->x_offset = ((I64)((point[0]/reproject_quantizer->x_scale_factor)/10000000))*10000000*reproject_quantizer->x_scale_factor;
        reproject_quantizer->y_offset = ((I64)((point[1]/reproject_quantizer->y_scale_factor)/10000000))*10000000*reproject_quantizer->y_scale_factor;
        reproject_quantizer->z_offset = ((I64)((point[2]/reproject_quantizer->z_scale_factor)/10000000))*10000000*reproject_quantizer->z_scale_factor;
      }

      int number_of_keys;
      GeoProjectionGeoKeys* geo_keys = 0;
      int num_geo_double_params;
      double* geo_double_params = 0;

      if (geoprojectionconverter.get_geo_keys_from_projection(number_of_keys, &geo_keys, num_geo_double_params, &geo_double_params, !geoprojectionconverter.has_projection(false)))
      {
        lasreader->header.set_geo_keys(number_of_keys, (LASvlr_key_entry*)geo_keys);
        free(geo_keys);
        if (geo_double_params)
        {
          lasreader->header.set_geo_double_params(num_geo_double_params, geo_double_params);
          free(geo_double_params);
        }
        else
        {
          lasreader->header.del_geo_double_params();
        }
        lasreader->header.del_geo_ascii_params();
      }
    }

    // do we need an extra pass

    BOOL extra_pass = laswriteopener.is_piped();

    // for piped output we need an extra pass

    if (extra_pass)
    {
      if (lasreadopener.is_piped())
      {
        fprintf(stderr, "ERROR: input and output cannot both be piped\n");
        usage(true);
      }


      if (verbose) fprintf(stderr, "extra pass for piped output: reading %lld points ...\n", lasreader->npoints);


      // maybe seek to start position

      if (subsequence_start) lasreader->seek(subsequence_start);

      while (lasreader->read_point())

      {
        if (lasreader->p_count > subsequence_stop) break;

        if (clip_to_bounding_box)
        {
          if (!lasreader->point.inside_box(lasreader->header.min_x, lasreader->header.min_y, lasreader->header.min_z, lasreader->header.max_x, lasreader->header.max_y, lasreader->header.max_z))
          {
            continue;
          }
        }

        if (reproject_quantizer)
        {
          lasreader->point.compute_coordinates();
          geoprojectionconverter.to_target(lasreader->point.coordinates);
          lasreader->point.compute_XYZ(reproject_quantizer);
        }
        lasinventory.add(&lasreader->point);
      }
      lasreader->close();

      lasreader->header.number_of_point_records = lasinventory.number_of_point_records;
      for (i = 0; i < 5; i++) lasreader->header.number_of_points_by_return[i] = lasinventory.number_of_points_by_return[i+1];
      if (reproject_quantizer) lasreader->header = *reproject_quantizer;
      lasreader->header.max_x = lasreader->header.get_x(lasinventory.max_X);
      lasreader->header.min_x = lasreader->header.get_x(lasinventory.min_X);
      lasreader->header.max_y = lasreader->header.get_y(lasinventory.max_Y);
      lasreader->header.min_y = lasreader->header.get_y(lasinventory.min_Y);
      lasreader->header.max_z = lasreader->header.get_z(lasinventory.max_Z);
      lasreader->header.min_z = lasreader->header.get_z(lasinventory.min_Z);

     // if (verbose) { fprintf(stderr,"extra pass took %g sec.\n", taketime()-start_time); start_time = taketime(); }

      if (verbose) fprintf(stderr, "piped output: reading %lld and writing %d points ...\n", lasreader->npoints, lasinventory.number_of_point_records);

    }
    else
    {
      if (reproject_quantizer)
      {
        saved_quantizer = new LASquantizer();
        *saved_quantizer = lasreader->header;
        lasreader->header = *reproject_quantizer;
      }

      //if (verbose) fprintf(stderr, "reading %lld and writing all surviving points ...\n", lasreader->npoints);

    }

    // check output

    if (!laswriteopener.active())
    {
      // create name from input name
      laswriteopener.make_file_name(lasreadopener.get_file_name());
    }

    // prepare the header for the surviving points

    strncpy(lasreader->header.system_identifier, "LAStools (c) by rapidlasso GmbH", 32);
    lasreader->header.system_identifier[31] = '\0';
    char temp[64];
    sprintf(temp, "las2las (version %d)", LAS_TOOLS_VERSION);
    strncpy(lasreader->header.generating_software, temp, 32);
    lasreader->header.generating_software[31] = '\0';


    LASwriter* laswriter = 0;
    // open laswriter
    if(is_mpi){
	// remove any existing out file, before opening with MPI_File_open
	if(rank==0){
	    remove(laswriteopener.get_file_name());
	}
	MPI_Barrier(MPI_COMM_WORLD);
    }


    laswriter = laswriteopener.open(&lasreader->header);
    if (laswriter == 0)
    {
         fprintf(stderr, "ERROR: could not open laswriter\n");
         byebye(true, argc==1);
    }
    // **************************************************************************************************
    if(is_mpi == 1){ // jdw, we do this because only rank 0 now writes the header in laswriter_las.cpp
      MPI_File fh = laswriter->get_MPI_File();
      MPI_Offset offset;
      //MPI_File_get_position(fh, &offset);
      //printf ("offset %lld, rank %i fh %lld\n", offset, rank, fh);
      if(rank==0){
           MPI_File_get_position(fh, &offset);
      }
      MPI_Bcast(&offset, 1, MPI_OFFSET, 0, MPI_COMM_WORLD);
      MPI_Barrier(MPI_COMM_WORLD);
      MPI_File_seek(fh, offset, MPI_SEEK_SET);

    }
    // ****************************************************************************************************



    // for piped output we need to re-open the input file

    if (extra_pass)
    {
      if (!lasreadopener.reopen(lasreader))
      {
        fprintf(stderr, "ERROR: could not re-open lasreader\n");
        byebye(true);
      }
    }
    else
    {
      if (reproject_quantizer)
      {
        lasreader->header = *saved_quantizer;
        delete saved_quantizer;
      }
    }

    // maybe seek to start position

    if (subsequence_start) lasreader->seek(subsequence_start);

    // loop over points

    if (point)
    {

      while (lasreader->read_point())

      {
        if (lasreader->p_count > subsequence_stop) break;

        if (clip_to_bounding_box)
        {
          if (!lasreader->point.inside_box(lasreader->header.min_x, lasreader->header.min_y, lasreader->header.min_z, lasreader->header.max_x, lasreader->header.max_y, lasreader->header.max_z))
          {
            continue;
          }
        }

        if (reproject_quantizer)
        {
          lasreader->point.compute_coordinates();
          geoprojectionconverter.to_target(lasreader->point.coordinates);
          lasreader->point.compute_XYZ(reproject_quantizer);
        }
        *point = lasreader->point;
        laswriter->write_point(point);
        // without extra pass we need inventory of surviving points
        if (!extra_pass) laswriter->update_inventory(point);
      }
      delete point;
      point = 0;
    }
    else // ***************************** MPI ********************************************************
    {
      // ***** Determine the start and stop points for this process *****
      I64 left_over_count = lasreader->npoints % process_count;
      I64 process_points = lasreader->npoints / process_count;
      subsequence_start = rank*process_points;
      subsequence_stop =  subsequence_start + process_points;
      if(rank == process_count-1) subsequence_stop += left_over_count;

      // ***** Set the input stream file offset for this process *****
      // subsequence_start parameter gets cast to U32 in the implementation of seek and overflows for large files
      // manually set the file offset instead for now
      //((LASreaderLAS*)lasreader)->stream->seek(subsequence_start);
      I64 header_end_read_position = lasreader->get_Stream()->tell();
      //printf("header end %lld subseqence_start * 28 %lld rank %i\n", header_end_read_position, subsequence_start*28, rank);
      lasreader->p_count = subsequence_start;
      lasreader->get_Stream()->seek(header_end_read_position + subsequence_start*28);
      //printf("seek pos first loop %lld rank %i\n", lasreader->get_Stream()->tell(), rank);


      if (verbose) fprintf(stderr, "reading %lli points, rank %i\n", subsequence_stop - subsequence_start, rank);

      // *****Read the file for the first time *****
      // this first read and filter of the file is to gather a count of points that pass the filter so that
      // write offsets can be set.
      I64 filtered_count = 0;
      while (lasreader->read_point()){
	  if (lasreader->p_count > subsequence_stop) break;
	  if (clip_to_bounding_box){
	    if (!lasreader->point.inside_box(lasreader->header.min_x, lasreader->header.min_y, lasreader->header.min_z, lasreader->header.max_x, lasreader->header.max_y, lasreader->header.max_z)){
	      continue;
	    }
	  }
	  filtered_count++;
      }
      // ***** Gather and set the write offset for this process *****
      I64* filtered_counts = (I64*)malloc(process_count * sizeof(I64));
      if(is_mpi)MPI_Barrier(MPI_COMM_WORLD);
      filtered_counts[rank] = filtered_count;
      if(is_mpi)MPI_Allgather(&filtered_count, 1, MPI_LONG_LONG, filtered_counts, 1, MPI_LONG_LONG, MPI_COMM_WORLD);
      if(is_mpi)MPI_Barrier(MPI_COMM_WORLD);

      if(debug) printf("filtered count %lli rank %i\n", filtered_counts[rank], rank);

      if(is_mpi)MPI_Barrier(MPI_COMM_WORLD);

      I64 write_point_offset = 0;
      for (int k=0; k < rank; k++){
	  write_point_offset += filtered_counts[k];
      }
      if(is_mpi){

        MPI_File fh = laswriter->get_MPI_File();
        MPI_Offset cur = 0;

        // jdw, todo, remove the hardcoding by adding methods to read point size from reader
        MPI_File_seek(fh, write_point_offset*28, MPI_SEEK_CUR);
        if(debug){
          MPI_File_get_position(fh, &cur);
          printf ("rank %i, write offset %lld\n", rank, write_point_offset*28);
        }
      }
      if(is_mpi)MPI_Barrier(MPI_COMM_WORLD);


      // ***** Read and filter the input file again, this time write the filtered point since output file offset in now known amd set *****
      //lasreader->seek(subsequence_start); // subsequence_start parameter gets cast to U32 in the implementation and overflows for large files
      // manually set the file offset instead for now
      //printf("header end %lld subseqence_start * 28 %lld rank %i\n", header_end_read_position, subsequence_start*28, rank);
      lasreader->p_count = subsequence_start;
      lasreader->get_Stream()->seek(header_end_read_position + subsequence_start*28);
      //printf("seek pos second loop %lld rank %i\n", lasreader->get_Stream()->tell(), rank);
      while (lasreader->read_point())
      {
          if (lasreader->p_count > subsequence_stop) break;

          if (clip_to_bounding_box)
          {
            if (!lasreader->point.inside_box(lasreader->header.min_x, lasreader->header.min_y, lasreader->header.min_z, lasreader->header.max_x, lasreader->header.max_y, lasreader->header.max_z))
            {
              continue;
            }
          }

          if (reproject_quantizer)
          {
            lasreader->point.compute_coordinates();
            geoprojectionconverter.to_target(lasreader->point.coordinates);
            lasreader->point.compute_XYZ(reproject_quantizer);
          }

          laswriter->write_point(&lasreader->point);
          // without extra pass we need inventory of surviving points
    	  if (!extra_pass){
            laswriter->update_inventory(&lasreader->point);
    	  }
      }
      //***** this is part of an mpi write optimization *****
      laswriter->get_Stream()->flushBytes();
    }

    // without the extra pass we need to fix the header now
    // ***** do the inventory reconciliation *****
    // ***** Reduce inventory information in rank 0 *****
    if (is_mpi){
        U32 number_of_point_records = 0;
        U32 number_of_points_by_return[8];
        for(int i = 0; i<8; i++)number_of_points_by_return[i] = 0;
        I32 max_X = 0;
        I32 min_X = 0;
        I32 max_Y = 0;
        I32 min_Y = 0;
        I32 max_Z = 0;
        I32 min_Z = 0;

        MPI_Reduce(&laswriter->inventory.number_of_point_records, &number_of_point_records, 1, MPI_UNSIGNED, MPI_SUM, 0, MPI_COMM_WORLD);
        MPI_Reduce(laswriter->inventory.number_of_points_by_return, number_of_points_by_return, 8, MPI_UNSIGNED, MPI_SUM, 0, MPI_COMM_WORLD);
        MPI_Reduce(&laswriter->inventory.max_X, &max_X, 1, MPI_INT, MPI_MAX, 0, MPI_COMM_WORLD);
        MPI_Reduce(&laswriter->inventory.min_X, &min_X, 1, MPI_INT, MPI_MIN, 0, MPI_COMM_WORLD);
        MPI_Reduce(&laswriter->inventory.max_Y, &max_Y, 1, MPI_INT, MPI_MAX, 0, MPI_COMM_WORLD);
        MPI_Reduce(&laswriter->inventory.min_Y, &min_Y, 1, MPI_INT, MPI_MIN, 0, MPI_COMM_WORLD);
        MPI_Reduce(&laswriter->inventory.max_Z, &max_Z, 1, MPI_INT, MPI_MAX, 0, MPI_COMM_WORLD);
        MPI_Reduce(&laswriter->inventory.min_Z, &min_Z, 1, MPI_INT, MPI_MIN, 0, MPI_COMM_WORLD);

        if (rank ==0){
            laswriter->inventory.number_of_point_records = number_of_point_records;
            for(int i=0; i<8; i++)laswriter->inventory.number_of_points_by_return[i] = number_of_points_by_return[i];
            laswriter->inventory.max_X = max_X;
            laswriter->inventory.min_X = min_X;
            laswriter->inventory.max_Y = max_Y;
            laswriter->inventory.min_Y = min_Y;
            laswriter->inventory.max_Z = max_Z;
            laswriter->inventory.min_Z = min_Z;
        }
    }

    if(rank == 0){
      if (!extra_pass)
      {
        if (reproject_quantizer) lasreader->header = *reproject_quantizer;
        laswriter->update_header(&lasreader->header, TRUE);
      }
    }
    if(is_mpi)MPI_Barrier(MPI_COMM_WORLD);
    if (verbose) { fprintf(stderr,"%lli surviving points written by rank: %i\n", laswriter->p_count, rank); }

    laswriter->close(FALSE);
    if(is_mpi)MPI_Barrier(MPI_COMM_WORLD);

    delete laswriter;
    lasreader->close();
    delete lasreader;
    if (reproject_quantizer) delete reproject_quantizer;

  }
  if(is_mpi)MPI_Finalize();

  time(&wall_end_time);

  if (verbose) { fprintf(stderr,"total time %.f sec, cpu time: %g sec. rank: %i\n", difftime(wall_end_time, wall_start_time), taketime()-start_time, rank); }
  return 0;
}
