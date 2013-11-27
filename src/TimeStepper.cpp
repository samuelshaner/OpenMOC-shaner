#include "TimeStepper.h"

/* TimeStepper constructor */
TimeStepper::TimeStepper(double start_time, double end_time, double dt_moc, double dt_cmfd){

    _start_time = start_time;
    _end_time = end_time;
    _dt_moc = dt_moc;
    _dt_cmfd = dt_cmfd;
    _state_times = new double[6];
    
    for (int i = 0; i < 6; i++)
	_state_times[i] = _start_time;
    
    _state_times[int(FORWARD)] += _dt_moc;
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


void TimeStepper::takeStep(){
    _state_times[int(PREVIOUS)] = _state_times[int(CURRENT)];
    _state_times[int(CURRENT)] += _dt_cmfd;
}


void TimeStepper::setTime(materialState state, double time){
    _state_times[(int)state] = time;
}


double TimeStepper::getImprovedRatio(){

    double ratio = (_state_times[int(CURRENT)] - _state_times [int(PREVIOUS_CONV)]) / _dt_moc;
    return ratio;
}

double TimeStepper::getTime(materialState state){
    return _state_times[(int)state];
}


double TimeStepper::getDtMOC(){
    return _dt_moc;
}


double TimeStepper::getDtCMFD(){
    return _dt_cmfd;
}


void TimeStepper::convergedMOCStep(){
    _state_times[int(PREVIOUS_CONV)] = _state_times[int(CURRENT)];
}

void TimeStepper::printTimes(){

    log_printf(NORMAL, "PREVIOUS_CONV TIME: %f", _state_times[int(PREVIOUS_CONV)]);
    log_printf(NORMAL, "PREVIOUS TIME: %f", _state_times[int(PREVIOUS)]);
    log_printf(NORMAL, "CURRENT TIME: %f", _state_times[int(CURRENT)]);
}


void TimeStepper::incrementTime(materialState state, double dt){
    _state_times[(int)state] += dt;
}
