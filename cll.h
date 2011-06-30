#ifndef ADVCL_CLL_H_INCLUDED
#define ADVCL_CLL_H_INCLUDED

#define __CL_ENABLE_EXCEPTIONS
#include "CL/cl.hpp"

// issue with using cl_float4 from cl_platform.h
// http://www.khronos.org/message_boards/viewtopic.php?f=28&t=1848
// typedef cl_float cl_float4 __attribute__ ((__vector_size__ (16), __may_alias__));
typedef struct Vec4
{
  float x,y,z,w;
    Vec4(){};
    //convenience functions
Vec4(float xx, float yy, float zz, float ww):
        x(xx),
        y(yy),
	z(zz),
        w(ww)
    {}
  void set(float xx, float yy, float zz, float ww) {
        x = xx;
        y = yy;
        z = zz;
	w = ww;
    }
} Vec4 __attribute__((aligned(16)));

typedef struct {
  Vec4 v;
  float f;
  float g;
  float h;
} Point;

class CL {
 public:
  //These are arrays we will use in this tutorial
  cl::Buffer cl_point;  
  // cl::Buffer cl_v;
  // cl::Buffer cl_fs;
  // cl::Buffer cl_gs;
  // cl::Buffer cl_hs;
  cl::Buffer cl_kf;  
  cl::Buffer cl_range;
  cl::Buffer cl_inc;
  cl::Buffer cl_xdim;
  cl::Buffer cl_ydim;
  cl::Buffer cl_zdim;
  cl::Buffer cl_binS;
  cl::Buffer cl_binsT;
  cl::Buffer cl_bins;
  cl::Buffer cl_binst;
  cl::Buffer cl_num;
  cl::Buffer cl_batchCount;
		
  int num;    //the number of particles
  size_t array_size; //the size of our arrays num * sizeof(Vec4)

  //default constructor initializes OpenCL context and automatically chooses platform and device
  CL();
  //default destructor releases OpenCL objects and frees device memory
  ~CL();

  //load an OpenCL program from a string
  void loadProgram(std::string kernel_source);
  //setup the data for the kernel

  void loadPoints(Point *, int, float *, float, int, int, int, int, int);

  void loadArguments();
  void pointKernel(float (*binsR)[110], float (*binstR)[110], float *, float *, int);
  void timing();

  // void loadData(Vec4 (*v)[8], float (*fs)[8],  float (*gs)[8], float (*hs)[8], int kappaFlag, float *, float, int, float (*bins)[110], float (*binst)[110]);
  // void loadPoints(Point *, int, float *, float, int, int, int, int, int);

  //void loadData(const std::vector<Vec4> &pos, float *);
  //these are implemented in part1.cpp (in the future we will make these more general)
  // void popCorn();
  //execute the kernel
  // void runKernel(float (*binsR)[110], float (*binstR)[110], float *, float *, int);

  //    real coders bare all
  //    private:

  unsigned int deviceUsed;
  std::vector<cl::Device> devices;
        
  cl::Context context;
  cl::CommandQueue queue;
  cl::Program program;
  cl::Kernel kernel;
  cl::Kernel initial;
  cl::Kernel summ;
        
  //debugging variables
  cl_int err;
  ///cl_event event;
  cl::Event event;
};

#endif
