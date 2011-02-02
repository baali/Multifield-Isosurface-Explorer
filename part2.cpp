#include <stdio.h>
#include "math.h"
#include "cll.h"
#include "util.h"

#include <vector>

// int flag = 0;
// int counter = 0;
// change 110 to proper values
void CL::loadData(Vec4 (*v)[8], float (*fs)[8], float (*gs)[8], float (*hs)[8], int kappaFlag, float range[],float inc, int numCells, float (*binS)[110], float (*binsT)[110])
{
    // for(int i = 0; i < 20; i++)
    //   printf("%f, %f, %f\n", v[0][i].x, v[0][i].y, v[0][i].z);
    num = numCells;
    //store the number of particles and the size in bytes of our arrays
    int v_size = num * 8 * sizeof(Vec4);
    int f_size = num * 8 * sizeof(float);
    int kf_size = sizeof(int);
    int range_size = sizeof(range);
    int inc_size = sizeof(float);
    int bins_size = num * 110 * sizeof(float);
    int binst_size = num * 110 * sizeof(float);

    try{
      cl_v = cl::Buffer(context, CL_MEM_READ_ONLY, v_size, NULL, &err);
      cl_fs = cl::Buffer(context, CL_MEM_READ_ONLY, f_size, NULL, &err);
      cl_gs = cl::Buffer(context, CL_MEM_READ_ONLY, f_size, NULL, &err);
      cl_hs = cl::Buffer(context, CL_MEM_READ_ONLY, f_size, NULL, &err);
      cl_kf = cl::Buffer(context, CL_MEM_READ_ONLY, kf_size, NULL, &err);
      cl_range = cl::Buffer(context, CL_MEM_READ_ONLY, range_size, NULL, &err);
      cl_inc = cl::Buffer(context, CL_MEM_READ_ONLY, inc_size, NULL, &err);
      cl_binS = cl::Buffer(context, CL_MEM_WRITE_ONLY, bins_size, NULL, &err);
      cl_binsT = cl::Buffer(context, CL_MEM_WRITE_ONLY, binst_size, NULL, &err);
      cl_bins = cl::Buffer(context, CL_MEM_WRITE_ONLY, 110 * sizeof(float), NULL, &err);
      cl_binst = cl::Buffer(context, CL_MEM_WRITE_ONLY, 110 * sizeof(float), NULL, &err);
      cl_num = cl::Buffer(context, CL_MEM_READ_ONLY, sizeof(int), NULL, &err);
    }
    catch (cl::Error er) {
        printf("ERROR: %s(%s)\n", er.what(), oclErrorString(er.err()));
    }

    // initialization of bins to zero
    initial = cl::Kernel(program, "initial", &err);

    try {
      err = queue.enqueueWriteBuffer(cl_binS, CL_TRUE, 0, bins_size, binS[0], NULL, &event);
      err = queue.enqueueWriteBuffer(cl_binsT, CL_TRUE, 0, binst_size, binsT[0], NULL, &event);
    }
    catch (cl::Error er) {
        printf("ERROR: %s(%s)\n", er.what(), oclErrorString(er.err()));
    }
    queue.finish();
    try {
      err = initial.setArg(0, cl_binS);
      err = initial.setArg(1, cl_binsT);
    }
    catch (cl::Error er) {
      printf("ERROR: %s(%s)\n", er.what(), oclErrorString(er.err()));
    }
    queue.finish();
    // kernel which does initialization of all zeros.
    try {    
      err = queue.enqueueNDRangeKernel(initial, cl::NullRange, cl::NDRange(num), cl::NullRange, NULL, &event); 
    }
    catch (cl::Error er) {
      printf("ERROR: %s(%s)\n", er.what(), oclErrorString(er.err()));
    }
    // printf("enqueueNDRangeKernel: %s\n", oclErrorString(err));
    queue.finish();

    // printf("Pushing data to the GPU\n");
    //push our CPU arrays to the GPU
    //data is tightly packed in std::vector starting with the adress of the first element
    try{
      err = queue.enqueueWriteBuffer(cl_v, CL_TRUE, 0, v_size, v[0], NULL, &event);
      err = queue.enqueueWriteBuffer(cl_fs, CL_TRUE, 0, f_size, fs[0], NULL, &event);
      err = queue.enqueueWriteBuffer(cl_gs, CL_TRUE, 0, f_size, gs[0], NULL, &event);
      err = queue.enqueueWriteBuffer(cl_hs, CL_TRUE, 0, f_size, hs[0], NULL, &event);
      err = queue.enqueueWriteBuffer(cl_kf, CL_TRUE, 0, kf_size, &kappaFlag, NULL, &event);
      err = queue.enqueueWriteBuffer(cl_range, CL_TRUE, 0, range_size, &range[0], NULL, &event);
      err = queue.enqueueWriteBuffer(cl_inc, CL_TRUE, 0, inc_size, &inc, NULL, &event);
    }
    catch (cl::Error er) {
        printf("ERROR: %s(%s)\n", er.what(), oclErrorString(er.err()));
    }
    queue.finish();
}

void CL::popCorn()
{
  // printf("in popCorn\n");
    //initialize our kernel from the program
    try{
      kernel = cl::Kernel(program, "part2", &err);
      summ = cl::Kernel(program, "summ", &err);
    }
    catch (cl::Error er) {
        printf("ERROR: %s(%s)\n", er.what(), oclErrorString(er.err()));
    }
    //set the arguements of our kernel
    try
    {
      err = kernel.setArg(0, cl_v);
      err = kernel.setArg(1, cl_fs);
      err = kernel.setArg(2, cl_gs);
      err = kernel.setArg(3, cl_hs);
      err = kernel.setArg(4, cl_kf);
      err = kernel.setArg(5, cl_range);
      err = kernel.setArg(6, cl_inc);
      err = kernel.setArg(7, cl_binS);
      err = kernel.setArg(8, cl_binsT);
    }
    catch (cl::Error er) {
      printf("ERROR: %s(%s)\n", er.what(), oclErrorString(er.err()));
    }
    
    //Wait for the command queue to finish these commands before proceeding
    queue.finish();
    // printf("Done popcorn\n");
}

void CL::runKernel(float (*binS)[110], float (*binsT)[110], float bins[110], float binst[110], int p)
{
    // printf("in runKernel\n");

    try {    
      err = queue.enqueueNDRangeKernel(kernel, cl::NullRange, cl::NDRange(p), cl::NullRange, NULL, &event); 
    }
    catch (cl::Error er) {
      printf("ERROR: %s(%s)\n", er.what(), oclErrorString(er.err()));
    }
    // printf("enqueueNDRangeKernel: %s\n", oclErrorString(err));
    queue.finish();
    // printf("done kernel\n");

    // err = queue.enqueueReadBuffer(cl_binS, CL_TRUE, 0, sizeof(float) * 110 * p, binS[0], NULL, &event);
    // err = queue.enqueueReadBuffer(cl_binsT, CL_TRUE, 0, sizeof(float) * 110 * p, binsT[0], NULL, &event);
    // printf("clEnqueueReadBuffer: %s\n", oclErrorString(err));
    // setting up summing kernel
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
    
    // for(int i = 0; i < p && flag == 0;i++)
    //   {
    // 	for(int j = 0;j < 100; j++)
    // 	  {
    // 	    if (isnan(binS[i][j]))
    // 	      {
    // 		flag = 1;
    // 		printf("nan at %d\n", j);
    // 		// break;
    // 	      }
    // 	    bins[j] += binS[i][j];
    // 	    binst[j] += binsT[i][j];
    // 	  }
    // 	if (flag)
    // 	  {
    // 	    printf("Cell number %d\n", counter+i);
    // 	    flag = 0;
    // 	  }	    
    //   }
    // counter += p;

    //clReleaseEvent(event);
}
