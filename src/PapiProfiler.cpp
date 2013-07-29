#include "PapiProfiler.h"

PapiProfiler::PapiProfiler(int num_threads, int num_codeSections) {

	if( num_threads > 0){
		_num_threads = num_threads;
		for(int i=0; i<num_threads; i++){
			_thrset_arr.push_back(papiThreadSet({PAPI_NULL, NULL, false}));
			_gbl_thrset_arr.push_back(papiThreadSet({PAPI_NULL, NULL, false}));
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
			_gbl_thrset_arr.push_back(papiThreadSet({PAPI_NULL, NULL, false}));
		}
	}
	if( num_threads < _num_threads){
		for(int i=0; i< (_num_threads - num_threads) ; i++){
			_thrset_arr.pop_back();
			_gbl_thrset_arr.pop_back();
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
		    	log_printf(ERROR,"Could not allocate thrset values");
		    }

	    	_gbl_thrset_arr[tid].values = 
		    	(long long *) malloc( sizeof(long long) * _EventCodes.size() );
		    if( _gbl_thrset_arr[tid].values == NULL){
		    	log_printf(ERROR,"Could not allocate gbl_thrset values");
		    }

		    for(i=0; i<_EventCodes.size(); i++){
		    	_thrset_arr[tid].values[i] = 0;
		    	_gbl_thrset_arr[tid].values[i] = 0;
		        retval = PAPI_add_event( *EventSet, _EventCodes[i]);
		        if ( retval != PAPI_OK )
		            handleError("Thread PAPI: Could not add event", retval);
		    }

		    _thrset_arr[tid].init = true;
		    _gbl_thrset_arr[tid].init = true;

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

void PapiProfiler::resetThreadCounts() {
	
	int i,j;
	for(i=0; i < _num_threads; i++){
		for(j=0; j < _EventCodes.size(); j++){
			_thrset_arr[i].values[j] = 0;
			_gbl_thrset_arr[i].values[j] = 0;
		}
	}

}

int PapiProfiler::deleteThreadSet() {
    
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

void PapiProfiler::gblAccum() {

	for(int i=0; i<_EventCodes.size(); i++){
		for(int j=0; j<_num_threads; j++){
			if(_gbl_thrset_arr[j].init == false ) 
				_gbl_thrset_arr[j].values[i] = 0;
			_gbl_thrset_arr[j].values[i] += _thrset_arr[j].values[i];
		}
	}

}

void PapiProfiler::printEventCountsPerUnit(int reduce, int perunit) {
    
    int i,j;
	int *accum = (int *) malloc(_EventCodes.size() * sizeof(int));
	double ratio;
    char eventname[CHARBUF];

    unsigned int n,d;

    log_printf(NORMAL,"Number of threads: %d", _num_threads);

    if( reduce == THR_SEPARATE){
		for(j=0; j<_EventCodes.size(); j++){
		    for(i=0; i<_num_threads; i++){
	    		PAPI_event_code_to_name(_EventCodes[j],eventname);
	    		log_printf(RESULT, "[%s]: TID [%d]: %f", 
	    			eventname, i, _thrset_arr[i].values[j]/(double)perunit);
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
    		log_printf(RESULT, "[%s]: Thread Sum: %f", 
	    			eventname, accum[j]/(double)perunit);
	    }
	}
	if( reduce == THR_SEPARATE_GBL ){
		for(j=0; j<_EventCodes.size(); j++){
		    for(i=0; i<_num_threads; i++){
	    		PAPI_event_code_to_name(_EventCodes[j],eventname);
	    		log_printf(RESULT, "GBL [%s]: TID [%d]: %f", 
	    			eventname, i, _gbl_thrset_arr[i].values[j]/(double)perunit);
	    	}
	    }
	}
	if( reduce == THR_REDUCE_GBL ){
		for(j=0; j<_EventCodes.size(); j++){
			accum[j] = 0;
		    for(i=0; i<_num_threads; i++){
		    	accum[j] += _gbl_thrset_arr[i].values[j];
	    	}
	    	PAPI_event_code_to_name(_EventCodes[j],eventname);
    		log_printf(RESULT, "GBL [%s]: Thread Sum: %f", 
    			eventname, (unsigned int)accum[j]/(double)perunit);
	    }

	}

	free(accum);	
}

void PapiProfiler::printEventCounts(int reduce) {
    
    int i,j;
	int *accum = (int *) malloc(_EventCodes.size() * sizeof(int));
	double ratio;
    char eventname[CHARBUF];

    unsigned int n,d;

    log_printf(NORMAL,"Number of threads: %d", _num_threads);

    if( reduce == THR_SEPARATE){
		for(j=0; j<_EventCodes.size(); j++){
		    for(i=0; i<_num_threads; i++){
	    		PAPI_event_code_to_name(_EventCodes[j],eventname);
	    		log_printf(RESULT, "[%s]: TID [%d]: %lld", 
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
    		log_printf(RESULT, "[%s]: Thread Sum: %lld", 
	    			eventname, accum[j]);
	    }
	}
	if( reduce == THR_SEPARATE_GBL ){
		for(j=0; j<_EventCodes.size(); j++){
		    for(i=0; i<_num_threads; i++){
	    		PAPI_event_code_to_name(_EventCodes[j],eventname);
	    		log_printf(RESULT, "GBL [%s]: TID [%d]: %lld", 
	    			eventname, i, _gbl_thrset_arr[i].values[j]);
	    	}
	    }
	}
	if( reduce == THR_REDUCE_GBL ){
		for(j=0; j<_EventCodes.size(); j++){
			accum[j] = 0;
		    for(i=0; i<_num_threads; i++){
		    	accum[j] += _gbl_thrset_arr[i].values[j];
	    	}
	    	PAPI_event_code_to_name(_EventCodes[j],eventname);
    		log_printf(RESULT, "GBL [%s]: Thread Sum: %lld", eventname, accum[j]);
	    }
	    // n = (unsigned int) accum[1];
	    // d = (unsigned int) accum[0];
	    // log_printf(NORMAL,"N/D: %f", (double)n/d);
	}

	free(accum);
}

long long PapiProfiler::getThreadEventCount(char* Event, int tid) {
	

	int index;
	int EventCode;
	int retval;

    if ( (retval = PAPI_event_name_to_code(Event, &EventCode)) != PAPI_OK){
        return handleError("Could not convert event name to code", retval);
    }
	index = getEventCodeIndex(EventCode);
	if( index != -1){
		return _gbl_thrset_arr[tid].values[index];
	} else {
		log_printf(NORMAL, "Event %s not found in current EventSet", Event);
	}
	return -1;
}

int PapiProfiler::getEventCodeIndex(int EventCode) {
	int index = -1;
	for(int i=0; i<_EventCodes.size(); i++) {
		if(_EventCodes[i] == EventCode )
			return i;
	}

	return index;
}