#include "PapiProfiler.h"

PapiProfiler::PapiProfiler(int num_threads, int num_codeSections) {

	if( num_threads > 0){
		_num_threads = num_threads;
		for(int i=0; i<num_threads; i++){
			_thrset_arr.push_back(papiThreadSet({PAPI_NULL, NULL, false}));
		}
	}
	else{
		log_printf(NORMAL, "Num threads must be positive integer");
	}

	_num_codeSections = num_codeSections;
	for(int i=0; i < num_codeSections; i++){
		_CodeSections.push_back(_EventCodes);
	}

	log_printf(NORMAL, "PAPI Profiler initialized");

}

PapiProfiler::~PapiProfiler() {

}

int PapiProfiler::init() {

	int retval;
	int status;

	status = PAPI_is_initialized();
	if( status == PAPI_NOT_INITED ){

	    /* Init PAPI library */
	    log_printf(NORMAL, "Intializing PAPI...");
	    if ((retval = PAPI_library_init(PAPI_VER_CURRENT)) != PAPI_VER_CURRENT) {
	        handleError("Could not init PAPI", retval);
	    }

	    /* Init threading */
	    retval = PAPI_thread_init(pthread_self);
	    if( retval != PAPI_OK ){
	        handleError("Could not init PAPI threads", retval);
	    }
	}

    return retval;

}

/**
 * @brief Adds an event to the Event Set
 * @details TODO: add detailed description
 * @param event the name of the event to add to the Event Set
 * @return PAPI error code
 */
int PapiProfiler::addEvent(char *event) {

    int EventCode;
    int retval;

    if ( (retval = PAPI_event_name_to_code(event, &EventCode)) != PAPI_OK){
        return handleError("Could not convert event name to code", retval);
    }
    _EventCodes.push_back(EventCode);

    for(int i=0; i<_CodeSections.size(); i++)
    	_CodeSections[i].push_back(EventCode);

    log_printf(NORMAL, "Added [%s]: %lu total events.", event, _EventCodes.size());
    
    return retval;
}

int PapiProfiler::getNumEvents() {
	return _EventCodes.size();
}

int PapiProfiler::clearEvents() {

	int retval;
    int num_events = _EventCodes.size();

    _EventCodes.clear();

    for(int i=0; i<_CodeSections.size(); i++)
    	_CodeSections[i].clear();

    log_printf(NORMAL, "Removed all %d events to be counted ", num_events);
    return retval;
}

void PapiProfiler::setNumThreads(int num_threads){

	if( num_threads > _num_threads){
		for(int i=0; i< (num_threads - _num_threads) ; i++){
			_thrset_arr.push_back(papiThreadSet({PAPI_NULL, NULL, false}));
		}
	}
	if( num_threads < _num_threads){
		for(int i=0; i< (_num_threads - num_threads) ; i++){
			_thrset_arr.pop_back();
		}
	}
	if( num_threads != _num_threads)
		log_printf(NORMAL,"Changed PAPI numthreads from %d to %d", 
			_num_threads, num_threads);
	_num_threads = num_threads;
}

int PapiProfiler::initThreadSet(int tid) {

    int i;
    int retval;
    int status;

    int *EventSet = &(_thrset_arr[tid].EventSet);

    if( _EventCodes.size() > 0 ){

	    if( *EventSet == PAPI_NULL ){
		    retval = PAPI_create_eventset(EventSet);
		    if ( retval != PAPI_OK ) {
		        handleError("Could not create Event Set for thread", retval);
		    }
		    log_printf(NORMAL, "Created EventSet %d", tid);
	    }

		if( _thrset_arr[tid].init == false ){

		    _thrset_arr[tid].values = 
		    	(long long *) malloc( sizeof(long long) * _EventCodes.size() );
		    if( _thrset_arr[tid].values == NULL){
		    	log_printf(ERROR,"Could not allocate thr set values");
		    }

		    for(i=0; i<_EventCodes.size(); i++){
		    	_thrset_arr[tid].values[i] = 0;
		        retval = PAPI_add_event( *EventSet, _EventCodes[i]);
		        if ( retval != PAPI_OK )
		            handleError("Thread PAPI: Could not add event", retval);
		    }

		    _thrset_arr[tid].init = true;

		    log_printf(NORMAL, "Papi Thread %d EventSet initialized", tid);

		}
	}

    return retval;
}

int PapiProfiler::startThreadSet(int tid) {

    int retval;
    int status;
    int *EventSet = &(_thrset_arr[tid].EventSet);

    if( _EventCodes.size() > 0 ){
	    if( (retval = PAPI_state(*EventSet, &status)) != PAPI_OK)
	    	handleError("Could not get EventSet status", retval);

	    if( status == PAPI_STOPPED ) {

		    retval = PAPI_start(*EventSet);
		    if(retval != PAPI_OK)
		        handleError("Thread PAPI: Could not start event", retval);
	        // log_printf(NORMAL, "PAPI Thread %d started counting", tid);
		}
	}

    return retval;
}

int PapiProfiler::accumThreadSet(int tid) {

	int retval;
	int status;
	int *EventSet = &(_thrset_arr[tid].EventSet);
	long long *values = _thrset_arr[tid].values;

	if( _EventCodes.size() > 0 ){
	    if( (retval = PAPI_state(*EventSet, &status)) != PAPI_OK)
	    	handleError("Could not get EventSet status", retval);

	    if( status == PAPI_RUNNING ) {
	    	retval = PAPI_accum(*EventSet, values);
		    if ( retval != PAPI_OK ) {
		        handleError("Thread PAPI: Could not accum events", retval);
		    }
	    }
	}	

	return retval;
}

int PapiProfiler::stopThreadSet(int tid) {

    int retval;
    int status;
    int *EventSet = &(_thrset_arr[tid].EventSet);
    long long *values = _thrset_arr[tid].values; 

    retval = PAPI_state(*EventSet, &status);

    if( _EventCodes.size() > 0 ){
	    if( status == PAPI_RUNNING ){
		    retval = PAPI_stop( *EventSet, NULL);
		    if ( retval != PAPI_OK ) {
		        handleError("Thread PAPI: Could not stop events", retval);
		    }
		} else {
			log_printf(ERROR, "EventSet %d not running", tid);
		}
	}

	// log_printf(NORMAL, "PAPI Thread %d stopped counting", tid);

    return retval;
}

int PapiProfiler::clearThreadSet() {
    
    int retval;
    int i;

    for(i=0; i < _num_threads; i++){
	    retval = PAPI_cleanup_eventset(_thrset_arr[i].EventSet);
	    if(retval != PAPI_OK) {
	        handleError("Thread PAPI: Could not clear event set", retval);
	    }
	    free(_thrset_arr[i].values);
	}

    return retval;
}

int PapiProfiler::handleError(const char *msg, int retval) {
    printf("%s: %s\n", msg, PAPI_strerror(retval));
    return retval;
}

void PapiProfiler::printEventCounts(int reduce) {
    
    int i,j;
	int *accum = (int *) malloc(_EventCodes.size() * sizeof(int));
	double l1_dcm_per_cycle;
    char eventname[CHARBUF];

    if( reduce == THR_SEPARATE){
		for(j=0; j<_EventCodes.size(); j++){
		    for(i=0; i<_num_threads; i++){
	    		PAPI_event_code_to_name(_EventCodes[j],eventname);
	    		log_printf(NORMAL, "[%s]: TID [%d]: %lld", 
	    			eventname, i, _thrset_arr[i].values[j]);
	    	}
	    }
	}
	if( reduce == THR_REDUCE ){
		for(j=0; j<_EventCodes.size(); j++){
			accum[j] = 0;
		    for(i=0; i<_num_threads; i++){
		    	accum[j] += _thrset_arr[i].values[j];
	    	}
	    	PAPI_event_code_to_name(_EventCodes[j],eventname);
    		log_printf(NORMAL, "[%s]: Thread Sum: %lld", 
	    			eventname, accum[j]);
	    }
	    l1_dcm_per_cycle = (double) accum[2]/ (double) accum[1];
	    log_printf(NORMAL, "[L1_DCM/CYCLE]: %f", l1_dcm_per_cycle);
	}

	free(accum);
}