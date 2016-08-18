p_las2las 

Background:

The LAStools las2las application and supporting LASlib library were extended 
with the MPI to allow the application to be run in parallel on a compute 
cluster. The goal was an application that would scale to arbitrarily large 
input, limited only by the amount of disk space needed to store the input and 
output file. No intermediate files are generated and individual process
memory needs are not determined by the size of the input or output.

The algorithm_and_test_results.pdf slide presentation and 
lidar_parallel_processing.pdf paper contain a description of the implementation
and test results on large LAS files run on the Stampede cluster at the 
University of Texas-Austin. 

Dependencies:
An MPI implementation must be installed. OpenMPI 1.6 and 1.8 are known to work
and were used in development and testing. mpic++ must be found in your PATH. 

Install:

git clone https://github.com/jwend/p_las2las
cd p_las2las
make

The p_las2las executable is in the bin directory.

Test:

mpirun -n 4 bin/p_las2las -i data/test.las -o out.las -keep_class 2 -verbose

Limitations and Supported Features:

p_las2las works only with LAS version 1.0, 1.1, and 1.2 and produces only 
corresponding LAS file output. Version 1.2 was tested most extensively with 
up to the 117 GB file size input. Most filters and re projection flags should
work with -keep_return, -keep_class, -keep_xy, -epsg, -target_longlat tested 
most extensively. For example here is a run command on a cluster with a 
111GB LAS 1.2 file over the Denver area:

mpirun -n 1024 ./p_las2las -i denver.las -o denver_clip_rtn1_r.las  -keep_xy 
510000 4396000 516000 4404000 -keep_return 1 -epsg 26913 -target_longlat 
-verbose







