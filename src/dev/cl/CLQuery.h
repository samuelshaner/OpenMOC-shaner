/**
 * @file CLQuery.h
 * @brief Routines to check machine for an NVIDIA GPU and print GPU 
 *        characteristics to the screen.
 * @author May 30, 2013
 * @author Robert Sloan, MIT, Course 8 (rsloan@mit.edu)
 */

#ifndef CLQUERY_H_
#define CLQUERY_H_

#include "CLInstance.hpp"
#include "../../log.h"

CLInstance * clQueueInst = (CLInstance*) 0;
bool machineContainsCL();

void printBasicCLInfo();
void printDetailedCLInfo();
int getNumThreadsInWorkgroup();

#endif /* CLQUERY_H_ */
