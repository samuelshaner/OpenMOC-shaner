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
    PREVIOUS_CONV,
    PREVIOUS,
    CURRENT,
    FORWARD,
    FSR,
    FSR_OLD,
    FORWARD_PREV
};

class TimeStepper {
private:
    double _start_time;
    double _end_time;
    double* _state_times;
    double _dt_moc;
    double _dt_cmfd;
    
public:
    TimeStepper(double start_time, double end_time, double dt_moc, double dt_cmfd);
    
    virtual ~TimeStepper();
    
    /* setters and workers */
    void takeStep();
    void setTime(materialState state, double time);
    void convergedMOCStep();
    void printTimes();
    void incrementTime(materialState state, double dt);
        
    /* getters */
    double getStartTime();
    double getEndTime();
    double getImprovedRatio();
    double getTime(materialState state);
    double getDtMOC();
    double getDtCMFD();
};

#endif
