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

///////////////////////////////////////////////////////////////
// CLArgument: Wrapper class for arguments to kernels
//

class CLArgument {
public:
    size_t argSize;
    void * arg;

    CLArgument(size_t _argSize, void* _arg);
    CLArgument(cl_int num);
    ~CLArgument();


    template<typename T>
    static CLArgument* toArg(T x) {
        T* stor = (T*)malloc(sizeof(T));
        *stor = x;
        return new CLArgument(sizeof(T), stor);
    }
};

///////////////////////////////////////////////////////////////
// Abbreviations for making arguments
//

typedef CLArgument CLArg;
typedef CLArgument CLA;

CLA* narg(cl_int num);
CLA* barg(size_t size, void* buf);

template<typename T>
CLA* targ(int length, T* buf) {
    return new CLA(length * sizeof(T), (void*)buf);
}

template<typename T>
CLArgument* toCLArg(T x) {
    return CLArgument::toArg<T>(x);
};

///////////////////////////////////////////////////////////////
// CLWorkDimensions: Class for specifying workload
//

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

// Abbreviations
typedef CLWorkDimensions CLDim;
typedef CLWorkDimensions CLWD;

///////////////////////////////////////////////////////////////
// CLInstance: Class for holding global instance
//

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

///////////////////////////////////////////////////////////////
// CLKernel: Wrapper class for OpenCL Kernels and execution
//

class CLKernel {
private:
    CLInstance parent;
    comp::program program;

public:
    comp::kernel kernel;
    std::string kernelName;

    CLKernel(CLInstance inst, comp::program _program, std::string name);
    ~CLKernel();

    void setKernelArgs(int count, ...);
    void setKernelArgVector(std::vector<CLArgument*> args);
    void setKernelArg(cl_uint argNum, CLArgument* arg);
};

///////////////////////////////////////////////////////////////
// CLMemory: Wrapper class for OpenCL Memory Buffers and
//           read/write operations
//

template<typename T>
class CLMemory {
private:
    CLInstance parent;
    comp::future<comp::buffer_iterator<T> > *
        iterator_waiting = 0;
    comp::future<T*> *
        buffer_waiting = 0;
    bool autosync;

public:
    comp::vector<T> vec;
    comp::buffer buf;

    CLMemory(CLInstance inst, size_t size, bool _autosync = true) :
        parent(inst),
        iterator_waiting(new comp::future<comp::buffer_iterator<T> >()),
        buffer_waiting(new comp::future<T *>()),
        autosync(_autosync),
        vec(comp::vector<T>(size, inst.ctx)),
        buf(vec.get_buffer())
    {};

    CLMemory(CLInstance inst, comp::vector<T> v, bool _autosync = true) :
        parent(inst),
        iterator_waiting(new comp::future<comp::buffer_iterator<T> >()),
        buffer_waiting(new comp::future<T *>()),
        autosync(_autosync),
        vec(v),
        buf(vec.get_buffer())
    {};

    ~CLMemory() {
        delete iterator_waiting;
        delete buffer_waiting;
    };


    void read(int hostBufLength, T * hostBuf) {
        *iterator_waiting =
            comp::copy_async(
                hostBuf,
                hostBuf + hostBufLength*sizeof(T),
                vec.begin(),
                parent.queue);

        if(autosync) wait();
    };

    void write(cl_uint hostBufLength, T * hostBuf) {
        *buffer_waiting =
            comp::copy_async(
                vec.begin(),
                vec.end(),
                hostBuf,
                parent.queue);

        if(autosync) wait();
    };

    void readVec(std::vector<T> hostVec) {
        // copy data from the host to the device
        *iterator_waiting =
            comp::copy_async(
                hostVec.begin(),
                hostVec.end(),
                vec.begin(),
                parent.queue);

        if(autosync) wait();
    };

    void writeVec(std::vector<T> hostVec) {
        // copy data from the device to the host
        *iterator_waiting =
            comp::copy_async(
                vec.begin(),
                vec.end(),
                hostVec.begin(),
                parent.queue);

        if(autosync) wait();
    };

    void wait() {
        buffer_waiting->wait();
        iterator_waiting->wait();
    };

    CLArg* arg() {
        cl_mem* clMemStor = (cl_mem*)malloc(sizeof(cl_mem));
        *clMemStor = buf.get();
        return new CLArgument(sizeof(cl_mem),clMemStor);
    };
};
