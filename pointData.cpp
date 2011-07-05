#include <stdio.h>
#include "math.h"
#include "cll.h"
#include "util.h"
#include <sys/time.h>

#include <vector>

timeval tim;
double t1, t2, t3;
double t_summ = 0;
double t_kappa = 0;
double t_transfer = 0;
int flag = 0;
int counter = 0;

void CL::timing()
{
  printf("Time spent in main loop %f\n", t_kappa);
  printf("Time spent in summation %f\n", t_summ);
  printf("Time spent in transfer %f\n", t_transfer);
}

int initFlag = 0;
// int pointCount = 1000000;
int step ;
void CL::loadPoints(Point *point, int kappaFlag, float range[],
		    float inc, int batchCount, int numCells,
		    int XDIM, int YDIM, int ZDIM)
{
  // This step would change based on how many blocks are left after last 
  // iteration.
  step = (numCells / ((XDIM - 1)*(YDIM - 1))) + 1;
  int initCellId = batchCount * (XDIM - 1) * (YDIM - 1) * (step - 1);
  int endCellId = initCellId + numCells - 1;
  // printf("Init cell id %d and end cell id %d\n", initCellId, endCellId);
  int initPointId = (initCellId % (XDIM - 1)) +
    ((initCellId / (XDIM - 1))%(YDIM - 1)) * (XDIM) +
    (initCellId / ((XDIM - 1) * (YDIM - 1))) * (XDIM) * (YDIM);
  int endPointId = (endCellId % (XDIM - 1)) +
    ((endCellId / (XDIM - 1))%(YDIM - 1)) * (XDIM) +
    (endCellId / ((XDIM - 1) * (YDIM - 1))) * (XDIM) * (YDIM);
  endPointId += (XDIM * YDIM) + XDIM ;
  int numPoints = endPointId - initPointId + 1;
  // printf("point ids are %d, %d, numPoints %d\n", initPointId, endPointId,
  // numPoints);
  num = numPoints;
  int point_size =  num * sizeof(Point);
  int kf_size = sizeof(int);
  int range_size = sizeof(range);
  int inc_size = sizeof(float);
  // printf("Size of bins created %d\n", sizeof(float) * 110 * numCells);
  int bins_size = numCells * 110 * sizeof(float);
  int binst_size = numCells * 110 * sizeof(float);

  if (initFlag == 0)
    {
      gettimeofday(&tim, NULL);
      t1=tim.tv_sec+(tim.tv_usec/1000000.0);
      try{
	cl_point = cl::Buffer(context, CL_MEM_READ_ONLY, point_size, NULL, &err);
	cl_kf = cl::Buffer(context, CL_MEM_READ_ONLY, kf_size, NULL, &err);
	cl_range = cl::Buffer(context, CL_MEM_READ_ONLY, range_size, NULL, &err);
	cl_inc = cl::Buffer(context, CL_MEM_READ_ONLY, inc_size, NULL, &err);
	cl_xdim = cl::Buffer(context, CL_MEM_READ_ONLY, sizeof(int), NULL, &err);
	cl_ydim = cl::Buffer(context, CL_MEM_READ_ONLY, sizeof(int), NULL, &err);
	cl_zdim = cl::Buffer(context, CL_MEM_READ_ONLY, sizeof(int), NULL, &err);
	cl_binS = cl::Buffer(context, CL_MEM_WRITE_ONLY, bins_size, NULL, &err);
	cl_binsT = cl::Buffer(context, CL_MEM_WRITE_ONLY, binst_size, NULL, &err);
	cl_bins = cl::Buffer(context, CL_MEM_WRITE_ONLY, 110 * sizeof(float), NULL, &err);
	cl_binst = cl::Buffer(context, CL_MEM_WRITE_ONLY, 110 * sizeof(float), NULL, &err);
	cl_num = cl::Buffer(context, CL_MEM_READ_ONLY, sizeof(int), NULL, &err);
	cl_batchCount = cl::Buffer(context, CL_MEM_READ_WRITE, sizeof(int), NULL, &err);
      }
      catch (cl::Error er) {
	printf("ERROR: %s(%s)\n", er.what(), oclErrorString(er.err()));
      }
      // gettimeofday(&tim, NULL);
      // t1=tim.tv_sec+(tim.tv_usec/1000000.0);
      try{
	// Again this step should change ...
	// Address of starting point id.
	err = queue.enqueueWriteBuffer(cl_point, CL_TRUE, 0, point_size, &point[initPointId], NULL, &event);
	// printf("Passed points\n");
	err = queue.enqueueWriteBuffer(cl_kf, CL_TRUE, 0, kf_size, &kappaFlag, NULL, &event);
	err = queue.enqueueWriteBuffer(cl_range, CL_TRUE, 0, range_size, &range[0], NULL, &event);
	err = queue.enqueueWriteBuffer(cl_inc, CL_TRUE, 0, inc_size, &inc, NULL, &event);
	err = queue.enqueueWriteBuffer(cl_xdim, CL_TRUE, 0, sizeof(int), &XDIM, NULL, &event);
	err = queue.enqueueWriteBuffer(cl_ydim, CL_TRUE, 0, sizeof(int), &YDIM, NULL, &event);
	err = queue.enqueueWriteBuffer(cl_zdim, CL_TRUE, 0, sizeof(int), &ZDIM, NULL, &event);
	err = queue.enqueueWriteBuffer(cl_batchCount, CL_TRUE, 0, sizeof(int), &batchCount, NULL, &event);
      }
      catch (cl::Error er) {
	printf("ERROR: %s(%s)\n", er.what(), oclErrorString(er.err()));
      }
      gettimeofday(&tim, NULL);
      t2 = tim.tv_sec+(tim.tv_usec/1000000.0);
      t_transfer += (t2 - t1);
    }
  else
    {
      // In second iteration, to pass just required data!
      try{
	err = queue.enqueueWriteBuffer(cl_point, CL_TRUE, 0, point_size, &point[0], NULL, &event);
      }
      catch (cl::Error er) {
	printf("ERROR: %s(%s)\n", er.what(), oclErrorString(er.err()));
      }
    }
  queue.finish();
}

void CL::loadArguments()
{
  // printf("Loading arguments\n");
  // initialize our kernel from the program
  try{
    kernel = cl::Kernel(program, "pointSet", &err);
    summ = cl::Kernel(program, "summ", &err);
  }
  catch (cl::Error er) {
    printf("ERROR: %s(%s)\n", er.what(), oclErrorString(er.err()));
  }
  // set the arguements of our kernel
  if (initFlag == 0)
    {
      gettimeofday(&tim, NULL);
      t1=tim.tv_sec+(tim.tv_usec/1000000.0);
      try
	{
	  err = kernel.setArg(0, cl_point);
	  err = kernel.setArg(1, cl_batchCount);
	  err = kernel.setArg(2, cl_kf);
	  err = kernel.setArg(3, cl_range);
	  err = kernel.setArg(4, cl_inc);
	  err = kernel.setArg(5, cl_xdim);
	  err = kernel.setArg(6, cl_ydim);
	  err = kernel.setArg(7, cl_zdim);
	  err = kernel.setArg(8, cl_binS);
	  err = kernel.setArg(9, cl_binsT);
	}
      catch (cl::Error er) {
	printf("ERROR: %s(%s)\n", er.what(), oclErrorString(er.err()));
      }
      gettimeofday(&tim, NULL);
      t2 = tim.tv_sec+(tim.tv_usec/1000000.0);
      t_transfer += (t2 - t1);
      // printf("Set arguments\n");
      // initFlag = 1;
    }
  else
    {
      try
  	{
  	  err = kernel.setArg(0, cl_point);
  	}
      catch (cl::Error er) {
  	printf("ERROR: %s(%s)\n", er.what(), oclErrorString(er.err()));
      }
    }
  //Wait for the command queue to finish these commands before proceeding
  queue.finish();
  // printf("Done loading arguments\n");
}

// void CL::pointKernel(float (*binS)[110], float (*binsT)[110], float bins[110], float binst[110], int p)
void CL::pointKernel(float bins[110], float binst[110], int p)
{
  // printf("Running kernel\n");
  gettimeofday(&tim, NULL);
  t1=tim.tv_sec+(tim.tv_usec/1000000.0);
  try {    
    err = queue.enqueueNDRangeKernel(kernel, cl::NullRange, cl::NDRange(p), cl::NullRange, NULL, &event); 
  }
  catch (cl::Error er) {
    printf("ERROR: %s(%s)\n", er.what(), oclErrorString(er.err()));
  }
  // printf("enqueueNDRangeKernel: %s\n", oclErrorString(err));
  queue.finish();
  gettimeofday(&tim, NULL);
  t2 = tim.tv_sec+(tim.tv_usec/1000000.0);
  t_kappa += (t2 - t1);

  gettimeofday(&tim, NULL);
  t1=tim.tv_sec+(tim.tv_usec/1000000.0);


  try {
    err = queue.enqueueWriteBuffer(cl_bins, CL_TRUE, 0, 110 * sizeof(float), &bins[0], NULL, &event);
    err = queue.enqueueWriteBuffer(cl_binst, CL_TRUE, 0, 110 * sizeof(float), &binst[0], NULL, &event);
    err = queue.enqueueWriteBuffer(cl_num, CL_TRUE, 0, sizeof(int), &p, NULL, &event);
  }
  catch (cl::Error er) {
      printf("ERROR: %s(%s)\n", er.what(), oclErrorString(er.err()));
  }

  try {
    err = summ.setArg(0, cl_binS);
    err = summ.setArg(1, cl_binsT);
    err = summ.setArg(2, cl_bins);
    err = summ.setArg(3, cl_binst);
    err = summ.setArg(4, cl_num);
  }
  catch (cl::Error er) {
    printf("ERROR: %s(%s)\n", er.what(), oclErrorString(er.err()));
  }
  queue.finish();
  try {    
  err = queue.enqueueNDRangeKernel(summ, cl::NullRange, cl::NDRange(110), cl::NullRange, NULL, &event); 
  }
  catch (cl::Error er) {
      printf("ERROR: %s(%s)\n", er.what(), oclErrorString(er.err()));
  }
  queue.finish();
  // printf("done summ\n");
  try {
    err = queue.enqueueReadBuffer(cl_bins, CL_TRUE, 0, sizeof(float) * 110 , &bins[0], NULL, &event);
    err = queue.enqueueReadBuffer(cl_binst, CL_TRUE, 0, sizeof(float) * 110 , &binst[0], NULL, &event);
  }
  catch (cl::Error er) {
      printf("ERROR: %s(%s)\n", er.what(), oclErrorString(er.err()));
  }
  queue.finish();

  // int id;
  // printf("Number of cells done %d\n", p);
  // try {    
  //   // err = queue.enqueueReadBuffer(cl_batchCount, CL_TRUE, 0, sizeof(int), &id, NULL, &event);
  //   printf("Size of bins read %d\n", sizeof(float) * 110 * p);
  //   err = queue.enqueueReadBuffer(cl_binS, CL_TRUE, 0, sizeof(float) * 110 * p, binS[0], NULL, &event);
  //   err = queue.enqueueReadBuffer(cl_binsT, CL_TRUE, 0, sizeof(float) * 110 *p, binsT[0], NULL, &event);
  // }
  // catch (cl::Error er) {
  //   printf("ERROR While reading: %s(%s)\n", er.what(), oclErrorString(er.err()));
  // }
  // queue.finish();
  // // printf("id is %d\n", id);
  // for (int i = 0; i < p; i++)
  //   {
  //     for (int j = 0; j < 110; j++)
  //       {
  //         bins[j] += binS[i][j];
  //         binst[j] += binsT[i][j];
  //       }
  //   }
  gettimeofday(&tim, NULL);
  t2 = tim.tv_sec+(tim.tv_usec/1000000.0);
  t_summ += (t2 - t1);

  //clReleaseEvent(event);
}
