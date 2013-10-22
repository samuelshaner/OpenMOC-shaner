/*
 * TimeStepper.h
 *
 *  Created on: September 16, 2013
 *      Author: samuelshaner
 */

#ifndef TIMESTEPPER_H_
#define TIMESTEPPER_H_

#define _USE_MATH_DEFINES
#include <math.h>
#include "log.h"

/**
 * Material states
 */
enum materialState {
	REFERENCE,
	PREVIOUS_CONV,
	PREVIOUS,
	CURRENT,
	FORWARD,
	ADJ
};

class TimeStepper {
private:
  double _start_time;
  double _end_time;
  double* _state_times;
  double _dt_outer;
  double _dt_intermediate;
  
public:
  TimeStepper(double start_time, double end_time, double dt_outer, double dt_intermediate);
  virtual ~TimeStepper();
  
  double getStartTime();
  double getEndTime();
  void takeOuterStep();
  void takeIntermediateStep();
  void convergedOuterStep();
  void setTime(materialState state, double time);
  double getImprovedRatio();
  double getTime(materialState state);
  double getDtOuter();
  double getDtIntermediate();
  void printTimes();
  void incrementTime(materialState state, double dt);
};

#endif
