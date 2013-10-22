#include "TimeStepper.h"

/* TimeStepper constructor */
TimeStepper::TimeStepper(double start_time, double end_time, double dt_outer, double dt_intermediate){

  _start_time = start_time;
  _end_time = end_time;
  _dt_outer = dt_outer;
  _dt_intermediate = dt_intermediate;
  _state_times = new double[6];
  
  for (int i = 0; i < 6; i++)
    _state_times[i] = _start_time;

  _state_times[4] += _dt_outer;
}


/* TimeStepper destructor */
TimeStepper::~TimeStepper(){
}


double TimeStepper::getStartTime(){
  return _start_time;
}


double TimeStepper::getEndTime(){
  return _end_time;
}


void TimeStepper::takeOuterStep(){
  return;
}


void TimeStepper::takeIntermediateStep(){
  _state_times[2] = _state_times[3];
}


void TimeStepper::convergedOuterStep(){
  _state_times[1] += _dt_outer;
  _state_times[2] = _state_times[1];
  _state_times[3] = _state_times[1];
  _state_times[4] = _state_times[1]+_dt_outer;
}

void TimeStepper::setTime(materialState state, double time){
  _state_times[(int)state] = time;
}


double TimeStepper::getImprovedRatio(){

  double ratio = (_state_times[3] - _state_times [1]) / _dt_outer;
  return ratio;
}

double TimeStepper::getTime(materialState state){
  return _state_times[(int)state];
}


double TimeStepper::getDtOuter(){
  return _dt_outer;
}


double TimeStepper::getDtIntermediate(){
  return _dt_intermediate;
}


void TimeStepper::printTimes(){

  log_printf(NORMAL, "REFERENCE TIME: %f", _state_times[0]);
  log_printf(NORMAL, "PREVIOUS_CONV TIME: %f", _state_times[1]);
  log_printf(NORMAL, "PREVIOUS TIME: %f", _state_times[2]);
  log_printf(NORMAL, "CURRENT TIME: %f", _state_times[3]);
  log_printf(NORMAL, "FORWARD TIME: %f", _state_times[4]);
  log_printf(NORMAL, "ADJOINT TIME: %f", _state_times[5]);
}


void TimeStepper::incrementTime(materialState state, double dt){
  _state_times[(int)state] += dt;
}
