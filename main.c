

/*********************************************************************************************

    This is public domain software that was developed by or for the U.S. Naval Oceanographic
    Office and/or the U.S. Army Corps of Engineers.

    This is a work of the U.S. Government. In accordance with 17 USC 105, copyright protection
    is not available for any work of the U.S. Government.

    Neither the United States Government, nor any employees of the United States Government,
    nor the author, makes any warranty, express or implied, without even the implied warranty
    of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE, or assumes any liability or
    responsibility for the accuracy, completeness, or usefulness of any information,
    apparatus, product, or process disclosed, or represents that its use would not infringe
    privately-owned rights. Reference herein to any specific commercial products, process,
    or service by trade name, trademark, manufacturer, or otherwise, does not necessarily
    constitute or imply its endorsement, recommendation, or favoring by the United States
    Government. The views and opinions of authors expressed herein do not necessarily state
    or reflect those of the United States Government, and shall not be used for advertising
    or product endorsement purposes.
*********************************************************************************************/


/****************************************  IMPORTANT NOTE  **********************************

    Comments in this file that start with / * ! or / / ! are being used by Doxygen to
    document the software.  Dashes in these comment blocks are used to create bullet lists.
    The lack of blank lines after a block of dash preceeded comments means that the next
    block of dash preceeded comments is a new, indented bullet list.  I've tried to keep the
    Doxygen formatting to a minimum but there are some other items (like <br> and <pre>)
    that need to be left alone.  If you see a comment that starts with / * ! or / / ! and
    there is something that looks a bit weird it is probably due to some arcane Doxygen
    syntax.  Be very careful modifying blocks of Doxygen comments.

*****************************************  IMPORTANT NOTE  **********************************/


#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <math.h>
#include <string.h>

#include "nvutility.h"

#include "shapefil.h"
#include "version.h"


/*

  Module Name:    build_swbd

  Programmer(s):  Jan C. Depner (PFM Software, area.based.editor AT gmail DOT com)

  Date Written:   July 07, 2013

  Purpose:        build_swbd reads the uncompressed, one-degree, Shuttle Radar Topography Mission (SRTM) Water Body Data
                  ESRI shape files and creates a simple compressed coastline file.  Data is stored in a compressed
                  (difference coded), bit-packed binary format.  The data is written as unsigned character arrays so
                  that there are no "endian" issues when reading the data on different architectures.  The final file is
                  in the following format:


                  Version:

                      128 characters ASCII


                  Header:

                      180 X 360 groups of three, 32 bit integers stored as character buffers.  Each group consists of the
                      address of the block of segments for the associated one-degree cell, the total number of segments in
                      the cell, and the total number of vertices in the cell.  The order of the 180 X 360 cells goes from
                      west to east, south to north beginning at -90/-180.  That is, the first cell is -90/-180, the next is 
                      -90/-179, etc. until we reach -90/179 at which point we go to -89/-180 and so on.  We add 90 to all
                      latitudes and 180 to all longitudes so that we can work in positive numbers so, in essence we really
                      go from 0/0 to 89/359.


                  Cell records:

                      Each cell consists of a number of segments.  Each segment in the cell consists of the following:


                          5 bits:        count bits
                          5 bits:        lon offset bits
                          5 bits:        lat offset bits
                          count bits:    count of vertices in the segment
                          18 bits:       lon bias + 2**17
                          18 bits:       lat bias + 2**17
                          26 bits:       start lon (times 100000)
                          25 bits:       start lat (times 100000)

                          count * (lon offset bits + lat offset bits):   lon and lat offsets (plus biases) from previous point


                      Some notes on the sizes - count bits is the number of bits needed to store the count of vertices in the
                      current segment.  This should never exceed 32.  If it does we've got a problem.  Lat and lon offset
                      bits is the number of bits needed to store the delta coded offsets between each lat and lon and the
                      previous lat and lon (or the start lat and lon) in the segment.  This also should never exceed 32.  The
                      lat and lon biases are values that we add to each lat and lon offset in order to make sure that we store
                      all offsets as positive values (I hate messing with sign extension ;-)  We use 18 bits to store these
                      because they should never exceed +-100000.  If they do then you have a line segment that is greater than
                      one degree and that is probably bogus.  We add 2**17 (131071) to this number before we store it in order
                      to ensure that it is always positive (stinkin' sign extension again ;-)  We store the start lat and lon
                      in 25 and 26 bits respectively because they are stored as positive integers times 100000 (range
                      0-17999999 and 0-35999999 respectively).  This gives us a resolution of about 1 meter at the equator.

                      In the olden days (when dinosaurs roamed the earth) I would have stored the start lat and lon as offsets
                      from the corner of the cell in order to save those few bits.  I would have also computed the number of
                      bits needed to store all of the bit counts.  This is what I did in 1981 with the original CIA WDBII
                      data and again in 1989 with the WVS data.  We used to have to worry about every bit back then.  Now I
                      find it much easier to understand if I keep that kind of logistical nightmare to a minimum ;-)  The
                      savings in storage are pretty minimal anyway.


  Caveats:        Requires shapelib version 1.2.10 or newer (many thanks to Frank Warmerdam for the library).

                  Some points are lost due to fuzziness in longitude in the upper latitudes.  Since we're trying to make a
                  coastline file this shouldn't be a big problem.  If you want to use the containers just use the
                  original SWBD shape files.  The output file from this program is about 327MB in size.  The original,
                  uncompressed shapefiles use up about 3.5GB.  Even compressed, the original files are 840MB in size.

                  The output file can be read using the read_coast function.  Note that this program and the read_coast
                  function use the variable type definitions in build_swbd_pfm_nvtypes.h.  That file was created in the early
                  1990's at the Naval Oceanographic Office and the type definitions are meant to be architecture independent.


  Arguments:      Input directory followed by output file name, for example:

                  build_swbd /data1/SWBDdata coast_swbd.ccl

*/


int32_t main (int32_t argc, char **argv)
{
  SHPHandle         shpHandle;
  SHPObject         *shape = NULL;
  FILE              *tfp, *fp, *ofp;
  int32_t           i, j, k, lnh, ln, lth, lt, ds, type, numShapes, numParts, total, diff_x[2], diff_y[2], num_vertices, segCount, *segx, *segy, file_ext;
  int32_t           input_file_count, percent, old_percent, address, offset, xoff, yoff, num_segments, range_x, range_y, count_bits, lon_offset_bits;
  int32_t           lat_offset_bits, size, bias_x, bias_y, pos, max_bias, lon_start, lon_end, lat_start, lat_end;
  uint8_t           start_segment = NVFalse, file_check = NVFalse, bad_flag = NVFalse;
  double            minBounds[4], maxBounds[4], lon, lat, cornerx[2], cornery[2], slon, slat;
  char              fname[512], dirname[512], shpname[512], version[128], outname[512], lathem, lonhem, dataset[6] = {'a', 'e', 'f', 'i', 'n', 's'};
  uint8_t           *buffer, head_buf[12];


  printf ("\n\n%s\n\n", VERSION);


  if (argc < 3)
    {
      fprintf (stderr, "Usage: %s INPUT_DIR OUTPUT_FILE\n", argv[0]);
      fprintf (stderr, "If the output file name does not have a .ccl extension it will be added.\n");
      exit (-1);
    }


  strcpy (dirname, argv[1]);


  /*  Initialize variables  */

  input_file_count = 0;
  total = 0;
  segx = NULL;
  segy = NULL;
  fp = NULL;
  lon = -999.0;
  lat = -999.0;


  /*  Make sure we don't have any old cell files hanging around in case we crashed previously.  */

  for (i = 0 ; i < 180 ; i++)
    {
      for (j = 0 ; j < 360 ; j++)
        {
          sprintf (fname, "cell_%03d_%03d", j, i);
          remove (fname);
        }
    }


  /*  Loop for both hemispheres.  */

  for (lnh = 0 ; lnh < 2 ; lnh++)
    {
      if (lnh)
        {
          lonhem = 'e';
          lon_start = 0;
          lon_end = 180;
        }
      else
        {
          lonhem = 'w';
          lon_start = 1;
          lon_end = 181;
        }


      /*  Loop throught the longitudes.  */

      for (ln = lon_start ; ln < lon_end ; ln++)
        {


          /*  Loop for both hemispheres.  */

          for (lth = 0 ; lth < 2 ; lth++)
            {
              if (lth)
                {
                  lathem = 'n';
                  lat_start = 0;
                  lat_end = 90;
                }
              else
                {
                  lathem = 's';
                  lat_start = 1;
                  lat_end = 91;
                }


              /*  Loop through the latitudes.  */

              for (lt = lat_start ; lt < lat_end ; lt++)
                {
                  file_check = NVFalse;


                  /*  Check to make sure we have a valid file before we open the output file.  */

                  for (ds = 0 ; ds < 6 ; ds++)
                    {
                      sprintf (shpname, "%s/%1c%03d%1c%02d%1c.shp", dirname, lonhem, ln, lathem, lt, dataset[ds]);


                      /*  Make sure the file exists before we try to open it with the shape library.  */

                      if ((tfp = fopen (shpname, "rb")) != NULL)
                        {
                          file_check = NVTrue;
                          file_ext = ds;
                          fclose (tfp);
                          break;
                        }
                    }


                  if (file_check)
                    {
                      /*  Figure out where the boundaries of the one degree cell are and build the output (temporary) filename.  */

                      if (lnh)
                        {
                          if (lth)
                            {
                              sprintf (fname, "cell_%03d_%03d", ln + 180, lt + 90);

                              cornerx[0] = (ln + 180) * 3600.0;
                              cornerx[1] = (ln + 181) * 3600.0;
                              cornery[0] = (lt + 90) * 3600.0;
                              cornery[1] = (lt + 91) * 3600.0;
                            }
                          else
                            {
                              sprintf (fname, "cell_%03d_%03d", ln + 180, -lt + 90);

                              cornerx[0] = (ln + 180) * 3600.0;
                              cornerx[1] = (ln + 181) * 3600.0;
                              cornery[0] = (-lt + 90) * 3600.0;
                              cornery[1] = (-lt + 91) * 3600.0;
                            }
                        }
                      else
                        {
                          if (lth)
                            {
                              sprintf (fname, "cell_%03d_%03d", -ln + 180, lt + 90);

                              cornerx[0] = (-ln + 180) * 3600.0;
                              cornerx[1] = (-ln + 181) * 3600.0;
                              cornery[0] = (lt + 90) * 3600.0;
                              cornery[1] = (lt + 91) * 3600.0;
                            }
                          else
                            {
                              sprintf (fname, "cell_%03d_%03d", -ln + 180, -lt + 90);

                              cornerx[0] = (-ln + 180) * 3600.0;
                              cornerx[1] = (-ln + 181) * 3600.0;
                              cornery[0] = (-lt + 90) * 3600.0;
                              cornery[1] = (-lt + 91) * 3600.0;
                            }
                        }


                      /*  Open the output file.  */

                      if ((fp = fopen (fname, "ab")) == NULL)
                        {
                          perror (fname);
                          exit (-1);
                        }


                      /*  Define the input shape file name.  */

                      sprintf (shpname, "%s/%1c%03d%1c%02d%1c.shp", dirname, lonhem, ln, lathem, lt, dataset[file_ext]);


                      input_file_count++;


                      /*  Initialize loop variables  */

                      segCount = 0;
                      percent = 0;
                      old_percent = -1;


                      /*  Open shape file  */

                      shpHandle = SHPOpen (shpname, "rb");

                      if (shpHandle == NULL)
                        {
                          perror (shpname);
                          exit (-1);
                        }


                      fprintf (stderr,"Reading %s                        \r", shpname);
                      fflush (stderr);


                      /*  Get shape file header info  */

                      SHPGetInfo (shpHandle, &numShapes, &type, minBounds, maxBounds);


                      /*  Read all shapes  */

                      bad_flag = NVFalse;
                      for (i = 0 ; i < numShapes ; i++)
                        {
                          shape = SHPReadObject (shpHandle, i);

                          total += shape->nVertices;


                          /*  Get all vertices  */

                          if (shape->nVertices >= 2)
                            {
                              for (j = 0, numParts = 1 ; j < shape->nVertices ; j++)
                                {
                                  start_segment = NVFalse;


                                  /*  Check for start of a new segment.  */

                                  if (!j && shape->nParts > 0) start_segment = NVTrue;


                                  /*  If the previous point was directly on a boundary it was probably a closure line (SWBD shape files
                                      are closed polygons that define areas of water) so we throw it out.  */

                                  if (bad_flag)
                                    {
                                      start_segment = NVTrue;
                                      bad_flag = NVFalse;
                                    }


                                  /*  Check for the start of a new segment inside a larger group of points (this would be a "Ring" point).  */

                                  if (numParts < shape->nParts && shape->panPartStart[numParts] == j)
                                    {
                                      start_segment = NVTrue;
                                      numParts++;
                                    }


                                  /*  Bias lat and lon by 90 and 180 so that all points are positive  */

                                  lon = shape->padfX[j] + 180.0;
                                  lat = shape->padfY[j] + 90.0;


                                  /*  Position in seconds to be compared with the cell boundaries.  */

                                  slon = lon * 3600.0;
                                  slat = lat * 3600.0;


                                  /*  Check for points (almost) exactly on any of the boundaries.  The longitudes get a bit fuzzy as we move
                                      farther away from the equator.  We may lose a point or two here or there but we're trying to make coastline
                                      not containers.  */

                                  if (fabs (slon - cornerx[0]) < 1.00000000000000015 || fabs (slon - cornerx[1]) < 1.00000000000000015 ||
                                      fabs (slat - cornery[0]) < 1.0 || fabs (slat - cornery[1]) < 1.0)
                                    {
                                      bad_flag = NVTrue;
                                    }
                                  else
                                    {
                                      /*  Damn boundary conditions!  */

                                      if (lon == 360.0) lon = 359.99999;


                                      /*  Start a new segment  */

                                      if (start_segment)
                                        {
                                          /*  Close last segment, start new segment  */

                                          if (segCount > 1)
                                            {
                                              fwrite (&segCount, sizeof (int32_t), 1, fp);

                                              for (k = 0 ; k < segCount ; k++)
                                                {
                                                  fwrite (&segx[k], sizeof (int32_t), 1, fp);
                                                  fwrite (&segy[k], sizeof (int32_t), 1, fp);
                                                }
                                            }

                                          segCount = 0;
                                        }


                                      /*  Allocate memory for the new point.  */

                                      segx = (int32_t *) realloc (segx, (segCount + 1) * sizeof (int32_t));
                                      if (segx == NULL)
                                        {
                                          perror ("Allocating segx memory");
                                          exit (-1);
                                        }

                                      segy = (int32_t *) realloc (segy, (segCount + 1) * sizeof (int32_t));
                                      if (segy == NULL)
                                        {
                                          perror ("Allocating segy memory");
                                          exit (-1);
                                        }


                                      /*  Add point to current segment  */

                                      segx[segCount] = NINT (lon * 100000.0);
                                      segy[segCount] = NINT (lat * 100000.0);


                                      /*  Increment the point counter.  */

                                      segCount++;
                                    }
                                }
                            }


                          /*  Destroy the shape object.  */

                          SHPDestroyObject (shape);
                        }


                      /*  Close out the last segment is it's not already closed.  */

                      if (segCount > 1)
                        {
                          fwrite (&segCount, sizeof (int32_t), 1, fp);

                          for (k = 0 ; k < segCount ; k++)
                            {
                              fwrite (&segx[k], sizeof (int32_t), 1, fp);
                              fwrite (&segy[k], sizeof (int32_t), 1, fp);
                            }

                          segCount = 0;
                        }


                      /*  Close the input file.  */

                      SHPClose (shpHandle);


                      /*  Close the temporary output file.  */

                      fclose (fp);
                    }
                }
            }
        }
    }


  /*  Free the segment memory.  */

  if (segx != NULL) free (segx);
  if (segy != NULL) free (segy);


  /*  Set the loop variables.  */

  percent = 0;
  old_percent = -1;
  total = 0;


  /*  Create the final output file name.  */

  strcpy (outname, argv[2]);
  if (strcmp (&outname[strlen (outname) - 4], ".ccl")) sprintf (outname, "%s.ccl", argv[2]);


  fprintf (stderr,"\n\n%s\n\n", outname);
  fflush (stderr);


  /*  Try to open the output file.  */

  if ((ofp = fopen (outname, "wb")) == NULL)
    {
      perror (outname);
      exit (-1);
    }


  /*  Write the header  */

  memset (version, 0, 128);
  sprintf (version, "%s\n", FILE_VERSION);
  fprintf(stderr,"%s\n",version);
  fflush (stderr);
  fwrite (version, 128, 1, ofp);


  /*  Initialize the header area  */

  for (i = 0 ; i < 180 ; i++)
    {
      for (j = 0 ; j < 360 ; j++)
        {
          offset = (i * 360 + j) * (3 * sizeof (int32_t)) + 128;

          address = 0;
          num_segments = 0;
          num_vertices = 0;

          fseek (ofp, offset, SEEK_SET);

          pos = 0;
          bit_pack (head_buf, pos, 8 * sizeof (int32_t), address); pos += (8 * sizeof (int32_t));
          bit_pack (head_buf, pos, 8 * sizeof (int32_t), num_segments); pos += (8 * sizeof (int32_t));
          bit_pack (head_buf, pos, 8 * sizeof (int32_t), num_vertices);

          fwrite (head_buf, 3 * sizeof (int32_t), 1, ofp);
        }
    }


  /*  Compute the maximum delta value.  */

  max_bias = (int32_t) (pow (2.0, 17.0) - 1.0);


  /*  Latitude loop.  */

  for (i = 0 ; i < 180 ; i++)
    {


      /*  Longitude loop.  */

      for (j = 0 ; j < 360 ; j++)
        {
          num_segments = 0;
          num_vertices = 0;


          /*  Define the input (temporary) file name.  */

          sprintf (fname, "cell_%03d_%03d", j, i);


          /*  Try to open the input file.  */

          if ((fp = fopen (fname, "rb")) != NULL)
            {

              /*  Compute the offset in the header at which to write the address, the number of segments, and the number of vertices.  */

              offset = (i * 360 + j) * (3 * sizeof (int32_t)) + 128;
              address = ftell (ofp);


              /*  Read the segment count from the input file.  */

              while (fread (&segCount, sizeof (int32_t), 1, fp))
                {
                  /*  Just in case we happened to write an empty (or single point) segment between files ;-)  */

                  if (segCount > 1)
                    {
                      num_vertices += segCount;
                      num_segments++;
                      total += segCount;


                      /*  Allocate memory for the segment.  */

                      segx = (int32_t *) calloc (segCount, sizeof (int32_t));
                      if (segx == NULL)
                        {
                          perror ("Allocating segx memory");
                          exit (-1);
                        }

                      segy = (int32_t *) calloc (segCount, sizeof (int32_t));
                      if (segy == NULL)
                        {
                          perror ("Allocating segy memory");
                          exit (-1);
                        }


                      /*  Compute the maximum difference between adjacent points in the segment.  */

                      diff_x[0] = 99999999;
                      diff_x[1] = -99999999;
                      diff_y[0] = 99999999;
                      diff_y[1] = -99999999;

                      for (k = 0 ; k < segCount ; k++)
                        {
                          if (!fread (&segx[k], sizeof (int32_t), 1, fp))
			    {
			      fprintf (stderr, "Bad return in file %s, function %s at line %d.  This should never happen!", __FILE__, __FUNCTION__, __LINE__ - 2);
			      fflush (stderr);
			      exit (-1);
			    }
			  
                          if (!fread (&segy[k], sizeof (int32_t), 1, fp))
			    {
			      fprintf (stderr, "Bad return in file %s, function %s at line %d.  This should never happen!", __FILE__, __FUNCTION__, __LINE__ - 2);
			      fflush (stderr);
			      exit (-1);
			    }

                          if (k)
                            {
                              diff_x[0] = MIN (segx[k] - segx[k - 1], diff_x[0]);
                              diff_x[1] = MAX (segx[k] - segx[k - 1], diff_x[1]);
                              diff_y[0] = MIN (segy[k] - segy[k - 1], diff_y[0]);
                              diff_y[1] = MAX (segy[k] - segy[k - 1], diff_y[1]);
                            }
                        }


                      bias_x = -diff_x[0];
                      bias_y = -diff_y[0];


                      if (bias_x > max_bias || bias_x < -max_bias)
                        {
                          fprintf (stderr, "\n\nlon bias out of range, terminating!\n\n");
                          fprintf (stderr, "%d %d %d\n", i,j,bias_x);
                          exit (-1);
                        }


                      if (bias_y > max_bias || bias_y < -max_bias)
                        {
                          fprintf (stderr, "\n\nlat bias out of range, terminating!\n\n");
                          fprintf (stderr, "%d %d %d\n", i,j,bias_y);
                          exit (-1);
                        }


                      range_x = diff_x[1] - diff_x[0];
                      range_y = diff_y[1] - diff_y[0];


                      if (!range_x) range_x = 1;
                      if (!range_y) range_y = 1;


                      /*  Compute the number of bits needed to store the data.  */

                      count_bits = int_log2 (segCount) + 1;
                      lon_offset_bits = int_log2 (range_x) + 1;
                      lat_offset_bits = int_log2 (range_y) + 1;


                      /*  Compute the size, in bytes, of the write buffer.  */

                      size = 5 + 5 + 5 + count_bits + lon_offset_bits + lat_offset_bits + 18 + 18 + 26 + 25 + 
                        (segCount - 1) * (lon_offset_bits + lat_offset_bits);

                      size = size / 8 + 1;


                      /*  Allocate the write buffer space.  */

                      buffer = (uint8_t *) calloc (1, size);

                      if (buffer == NULL)
                        {
                          perror ("Allocating buffer");
                          exit (-1);
                        }


                      /*  Bit pack the data into the write buffer.  */

                      pos = 0;
                      bit_pack (buffer, pos, 5, count_bits); pos += 5;
                      bit_pack (buffer, pos, 5, lon_offset_bits); pos += 5;
                      bit_pack (buffer, pos, 5, lat_offset_bits); pos +=5;
                      bit_pack (buffer, pos, count_bits, segCount); pos += count_bits;
                      bit_pack (buffer, pos, 18, bias_x + max_bias); pos += 18;
                      bit_pack (buffer, pos, 18, bias_y + max_bias); pos += 18;
                      bit_pack (buffer, pos, 26, segx[0]); pos += 26;
                      bit_pack (buffer, pos, 25, segy[0]); pos += 25;


                      for (k = 1 ; k < segCount ; k++)
                        {
                          xoff = (segx[k] - segx[k - 1]) + bias_x;
                          yoff = (segy[k] - segy[k - 1]) + bias_y;

                          bit_pack (buffer, pos, lon_offset_bits, xoff); pos += lon_offset_bits;
                          bit_pack (buffer, pos, lat_offset_bits, yoff); pos += lat_offset_bits;
                        }


                      /*  Now, write the buffer to the output file.  */

                      fwrite (buffer, size, 1, ofp);


                      /*  Free the buffer and segment memory.  */

                      free (buffer);
                      free (segx);
                      free (segy);
                    }
                }


              /*  Close the input file and delete it.  */

              fclose (fp);
              remove (fname);


              /*  Write the address, number of segments, and number of vertices in the header  */

              fseek (ofp, offset, SEEK_SET);

              k = 8 * sizeof (int32_t);

              pos = 0;
              bit_pack (head_buf, pos, k, address); pos += k;
              bit_pack (head_buf, pos, k, num_segments); pos += k;
              bit_pack (head_buf, pos, k, num_vertices);

              fwrite (head_buf, 3 * sizeof (int32_t), 1, ofp);

              fseek (ofp, 0, SEEK_END);
            }
        }

      percent = (int32_t) (((float) i / 181.0) * 100.0);
      if (percent != old_percent)
        {
          fprintf (stderr, "%03d%% packed\r", percent);
          fflush (stderr);
          old_percent = percent;
        }
    }


  /*  Close the output file.  */

  fclose (ofp);


  fprintf (stderr, "100%% packed\n\n");
  fprintf (stderr, "Total points packed = %d\n\n", total);
  fflush (stderr);

  return (0);
}
