#include <stdio.h>
#include <string>
#include <iostream>

#include "cll.h"
#include "util.h"

CL::CL()
{
    // printf("Initialize OpenCL object and context\n");
    //setup devices and context
    std::vector<cl::Platform> platforms;
    err = cl::Platform::get(&platforms);
    // printf("cl::Platform::get(): %s\n", oclErrorString(err));
    // printf("platforms.size(): %d\n", platforms.size());

    deviceUsed = 0;
    err = platforms[0].getDevices(CL_DEVICE_TYPE_GPU, &devices);
    // printf("getDevices: %s\n", oclErrorString(err));
    // printf("devices.size(): %d\n", devices.size());
    int t = devices.front().getInfo<CL_DEVICE_TYPE>();
    // printf("type: device: %d CL_DEVICE_TYPE_GPU: %d \n", t, CL_DEVICE_TYPE_GPU);
    cl_context_properties properties[] = 
        { CL_CONTEXT_PLATFORM, (cl_context_properties)(platforms[0])(), 0};

    context = cl::Context(CL_DEVICE_TYPE_GPU, properties);
    devices = context.getInfo<CL_CONTEXT_DEVICES>();

    //create the command queue we will use to execute OpenCL commands
    try{
        queue = cl::CommandQueue(context, devices[deviceUsed], 0, &err);
    }
    catch (cl::Error er) {
        printf("ERROR: %s(%d)\n", er.what(), er.err());
    }

}

CL::~CL()
{
}


void CL::loadProgram(std::string source_str)
{
    // Program Setup
    //size_t program_length;
    // printf("load the program\n");
    // printf("%s\n", source_str.c_str());
    // printf("kernel size: %d\n", source_str.length()+1);
    // printf("kernel: \n %s\n", kernel_source.c_str());
    try
    {
        cl::Program::Sources source(1,
				    std::make_pair(source_str.c_str(), source_str.length()+1));
        program = cl::Program(context, source);
    }
    catch (cl::Error er) {
        printf("ERROR: %s(%s)\n", er.what(), oclErrorString(er.err()));
    }

    // printf("building program...\n");
    int t = devices.front().getInfo<CL_DEVICE_TYPE>();

    try
    {
      err = program.build(devices);
    }
    catch (cl::Error er) {
        printf("program.build: %s\n", oclErrorString(er.err()));
    //if(err != CL_SUCCESS){
    }
    catch (...) {
      printf("Some other error\n");
    }

    // printf("done building program\n");
    // std::cout << "Build Status: " << program.getBuildInfo<CL_PROGRAM_BUILD_STATUS>(devices[0]) << std::endl;
    // std::cout << "Build Options:\t" << program.getBuildInfo<CL_PROGRAM_BUILD_OPTIONS>(devices[0]) << std::endl;
    // std::cout << "Build Log:\t " << program.getBuildInfo<CL_PROGRAM_BUILD_LOG>(devices[0]) << std::endl;
}

