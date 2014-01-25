#include "CLQuery.h"

/**
 * @brief Queries a node to determine whether it contains one or more CLs.
 * @return True if the node contains a CL, false otherwise.
 */
bool machineContainsCL() {
    if(!clQueueInst) {
        clQueueInst = new CLInstance()
    }

    return (clQueueInst.devices.size() > 0);
}


/**
 * @brief Prints the basic device info for the CUDA-enabled device with ID=0.
 * @details Prints the name, compute capability, # multiprocessors and
 *          the clock rate of the device.
 */
void printBasicCLInfo() {

    if (!machineContainsCL()) {
        log_printf(WARNING, "Unable to print basic device info since no"
		 "opencl device is attached to the machine");
        return;
    }

    log_printf(NORMAL, "Device name: %s", clQueueInst.currentDevice.name());
    log_printf(NORMAL, "Device vendor: %s", clQueueInst.currentDevice.vendor());
    log_printf(NORMAL, "Device compute units: %d",
        clQueueInst.currentDevice.compute_units());
}


/**
 * @brief Prints the detailed device info for the CUDA-enabled device with ID=0.
 * @details Prints the total global and constant memory, shared memory and 
 *          registers per multiprocessor, # threads per warp, maximum # 
 *          threads per multiprocessor, maximum # threads per block, maximum
 *          threadblock dimensions, and maximum grid dimensions.
 */
void printDetailedCLInfo() {
    if (!machineContainsCL()) {
        log_printf(WARNING, "Unable to print basic device info since no"
		 "opencl device is attached to the machine");
        return;
    }

    log_printf(NORMAL, "Device name: %s", clQueueInst.currentDevice.name());
    log_printf(NORMAL, "Device vendor: %s", clQueueInst.currentDevice.vendor());
    log_printf(NORMAL, "Device compute units: %d",
        clQueueInst.currentDevice.compute_units());
    log_printf(NORMAL, "Device clock frequency: %d", 
	   clQueueInst.currentDevice.clock_frequency());
    log_printf(NORMAL, "Device local memory: %d", 
       clQueueInst.currentDevice.local_memory_size());
    log_printf(NORMAL, "Device global memory: %d", 
       clQueueInst.currentDevice.global_memory_size());
    log_printf(NORMAL, "Device max work group: %d", 
       clQueueInst.currentDevice.max_work_group_size());
}



/**
 * @brief Returns the number of threads in a workgroup for the attached CLInstance.
 */
int getNumThreadsInWorkgroup() {
    if (!machineContainsCL()) {
        log_printf(WARNING, "Unable to return the number of threads per workgroup "
		   "since no opencl device is attached to the machine");
        return 0;
    }

    return (int)clQueueInst.currentDevice.max_work_group_size()
}
