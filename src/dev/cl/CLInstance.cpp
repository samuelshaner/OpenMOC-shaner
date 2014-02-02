#include "CLInstance.hpp"

CLKernel::CLKernel(CLInstance inst, comp::program _program, std::string name) {
	parent = inst;
	kernelName = name;
	program = _program;

    kernel = program.create_kernel(name.c_str());
}

CLKernel::~CLKernel() {
}

void
CLKernel::setKernelArgs(int count, ...) {
	va_list arglist;
    std::vector<CLArgument*> kArgs;

    va_start(arglist, count);
    for(int j = 0; j < count; j++) {
        kArgs.push_back(va_arg(arglist, CLArgument*));
    }

    va_end(arglist);
    setKernelArgVector(kArgs);
}

void
CLKernel::setKernelArgVector(std::vector<CLArgument*> args) {
	for(cl_uint i=0; i < (cl_uint)args.size(); i++) {
		setKernelArg(i, args[i]);
    }
}

void
CLKernel::setKernelArg(cl_uint argNum, CLArgument* arg) {
    kernel.set_arg(argNum, arg->argSize, arg->arg);
}

CLArgument::CLArgument(size_t _argSize, void* _arg) {
	argSize = _argSize;
	arg = _arg;
}

CLArgument::CLArgument(cl_int num) {
    cl_int* nstor = (cl_int*)malloc(sizeof(cl_int));
    *nstor = num;
    CLArgument(sizeof(cl_int), nstor);
}

CLArgument::~CLArgument() {
    free(arg);
}

CLA* narg(cl_int num) {
    return new CLA(num);
}

CLA* barg(size_t size, void* buf) {
    return new CLA(size, buf);
}

CLWorkDimensions::CLWorkDimensions(
	size_t xLocal, size_t xGlobal) {

	size_t globalWork[3] = {xGlobal, 0, 0};
    size_t localWork[3] = {xLocal, 0, 0};

    numDims = 1;

    globalDims = (size_t*)globalWork;
    localDims = (size_t*)localWork;
}

CLWorkDimensions::CLWorkDimensions(
    size_t xLocal, size_t yLocal,
    size_t xGlobal, size_t yGlobal) {

	size_t globalWork[3] = {xGlobal, yGlobal, 0};
    size_t localWork[3] = {xLocal, yLocal, 0};

    numDims = 2;

    globalDims = (size_t*)globalWork;
    localDims = (size_t*)localWork;
}

CLWorkDimensions::CLWorkDimensions(
    size_t xLocal, size_t yLocal, size_t zLocal,
    size_t xGlobal, size_t yGlobal, size_t zGlobal) {

	size_t globalWork[3] = {xGlobal, yGlobal, zGlobal};
    size_t localWork[3] = {xLocal, yLocal, zLocal};

    numDims = 3;

    globalDims = (size_t*)globalWork;
    localDims = (size_t*)localWork;
}

CLInstance::CLInstance() :
	ctx(comp::system::default_context()),
	queue(comp::command_queue(
		ctx, comp::system::default_device(),
		comp::command_queue::enable_profiling)),
    platforms(comp::system::platforms()),
    devices(comp::system::devices()),
    currentDevice(comp::system::default_device())
{}

CLInstance::~CLInstance() { }

void
CLInstance::setCurrentDevice(int deviceIndex) {
	currentDevice = devices[deviceIndex];
}

void
CLInstance::wait() {
	queue.finish();
}

void
CLInstance::runOnce(CLKernel ker) {
	queue.enqueue_task(ker.kernel);
}

void
CLInstance::runOnceWithArgs(CLKernel ker, int count, ...) {
	va_list arglist;
    std::vector<CLArgument*> kArgs;

    va_start(arglist, count);
    for(int j = 0; j < count; j++) {
        kArgs.push_back(va_arg(arglist, CLArgument*));
    }

    va_end(arglist);
	runOnceWithArgVec(ker, kArgs);	
}

void
CLInstance::runOnceWithArgVec(CLKernel ker, std::vector<CLArgument*> args) {
	ker.setKernelArgVector(args);	
	queue.enqueue_task(ker.kernel);
}

void
CLInstance::run(CLKernel ker, CLWorkDimensions wd, int count, ...) {
	va_list arglist;
    std::vector<CLArgument*> kArgs;

    va_start(arglist, count);
    for(int j = 0; j < count; j++)
        kArgs.push_back(va_arg(arglist, CLArgument*));

    va_end(arglist);
    runWithArgVec(ker, wd, kArgs);
}

void
CLInstance::runWithArgVec(
	CLKernel ker,
	CLWorkDimensions wd,
	std::vector<CLArgument*> args) {

	ker.setKernelArgVector(args);

	queue.enqueue_nd_range_kernel(
		ker.kernel,
		wd.numDims,
		wd.workOffset,
		wd.globalDims,
		wd.localDims);
}

comp::program
CLInstance::createProgOnline(std::string fname) {
    return comp::program::create_with_source_file(fname, ctx);
}

comp::program
CLInstance::createProgOffline(std::string fname) {
    return comp::program::create_with_binary_file(fname, ctx);
}
