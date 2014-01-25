#pragma once 

#include <stdlib.h>
#include <stdarg.h>
#include <string>
#include <vector>
#include <map>

#include <CL/cl.h>

#include <boost/compute.hpp>

namespace comp = boost::compute;

typedef comp::program CLProgram;

class CLArgument {
public:
    size_t argSize;
    void * arg;

    CLArgument(size_t _argSize, void* _arg);
    ~CLArgument();
};

class CLWorkDimensions {
public:
    cl_uint numDims;
    size_t * localDims;
    size_t * globalDims;
    size_t workOffset[3] = {0,0,0};

    CLWorkDimensions(size_t xLocal, size_t xGlobal);
    CLWorkDimensions(
        size_t xLocal, size_t yLocal,
        size_t xGlobal, size_t yGlobal);
    CLWorkDimensions(
        size_t xLocal, size_t yLocal, size_t zLocal,
        size_t xGlobal, size_t yGlobal, size_t zGlobal);
};

// Forward Declaration
class CLKernel;

class CLInstance {
public:
    CLInstance();
    ~CLInstance();

    comp::context ctx;
    comp::command_queue queue;

    std::vector<comp::platform> platforms;
    std::vector<comp::device> devices;

    comp::device currentDevice;
    void setCurrentDevice(int deviceIndex);

    bool isError(cl_int retCode);
    void wait();

    void runOnce(CLKernel ker);
    void runOnceWithArgs(CLKernel ker, int count, ...);
    void runOnceWithArgVec(CLKernel ker, std::vector<CLArgument*> args); 
    void run(CLKernel ker, CLWorkDimensions wd, int count, ...);
    void runWithArgVec(CLKernel ker, CLWorkDimensions wd, std::vector<CLArgument*> args);

    CLProgram createProgOnline(std::string fname);
    CLProgram createProgOffline(std::string fname);
};

class CLKernel {
private:
    CLInstance parent;
    comp::program program;

public:
    comp::kernel kernel;
    std::string kernelName;

    CLKernel(CLInstance inst, comp::program _program, std::string name);
    ~CLKernel() {};

    void setKernelArgs(int count, ...);
    void setKernelArgVector(std::vector<CLArgument*> args);
    void setKernelArg(cl_uint argNum, CLArgument* arg);
};

class CLMemory {
private:
    CLInstance parent;

public:
    comp::buffer buf;

    CLMemory(CLInstance inst, size_t size);
    ~CLMemory() {};

    CLArgument* getAsArg();

    cl_int read(cl_uint hostBufSize, void * hostBuf);
    cl_int write(cl_uint hostBufSize, void * hostBuf);
};

template<typename T>
CLArgument* toArg(T x) {
    return new CLArgument(sizeof(T), &x);
}
