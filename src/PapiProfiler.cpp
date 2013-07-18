#include "PapiProfiler.h"

PapiProfiler::PapiProfiler(int num_threads) {

	if( num_threads > 0){
		_thrset_arr = (papiThreadSet *) malloc( sizeof(papiThreadSet) * num_threads);
		if( _thrset_arr == NULL){
			log_printf(ERROR, "Could not allocate thrset_arr");
			exit(1);
		}
		_num_threads = num_threads;
	}
	else{
		log_printf(ERROR, "Num threads must be positive integer");
	}

	log_printf(NORMAL, "PAPI Profiler created");

}

PapiProfiler::~PapiProfiler() {

	free(_thrset_arr);

}

int PapiProfiler::init() {

	int retval;

    /* Init PAPI library */
    printf("\nIntializing PAPI...\n");
    if ((retval = PAPI_library_init(PAPI_VER_CURRENT)) != PAPI_VER_CURRENT) {
        handleError("Could not init PAPI", retval);
    }

    retval = PAPI_thread_init(pthread_self);
    if( retval != PAPI_OK ){
        handleError("Could not init PAPI threads", retval);
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
    printf("Added [%s] to the Event Set. %lu total events.\n", event, _EventCodes.size());
    
    return retval;
}

int PapiProfiler::clearEvents() {

	int retval;

    int num_events = _EventCodes.size();
    _EventCodes.clear();
    printf("Removed %d events.\n", num_events);
    return retval;
}

void PapiProfiler::setNumThreads(int num_threads){

	if( num_threads > _num_threads){
		_thrset_arr = 
			(papiThreadSet *) realloc( _thrset_arr, sizeof(papiThreadSet) * num_threads);
		if( _thrset_arr == NULL){
			log_printf(ERROR,"Could not realloc thrset_arr");
		}

	}
	log_printf(NORMAL,"Changing PAPI numthreads from %d to %d", _num_threads, num_threads);
	_num_threads = num_threads;
}

/* Each thread calls this in the parallel section */
int PapiProfiler::initThreadSet() {

    int i,j;
    int retval;
    std::vector<int>::iterator it;

    for(j=0; j<_num_threads; j++){

	    int *EventSet = &(_thrset_arr[j].EventSet);
	    *EventSet = PAPI_NULL;

	    retval = PAPI_create_eventset(EventSet);
	    if ( retval != PAPI_OK ) {
	        handleError("Could not create Event Set for thread", retval);
	    }

	    for(i=0; i<_EventCodes.size(); i++){
	        retval = PAPI_add_event( *EventSet, _EventCodes[i]);
	        if ( retval != PAPI_OK ) {
	            handleError("Thread PAPI: Could not add event", retval);
	        }
	    }

	    _thrset_arr[j].values = 
	    	(long long *) malloc( sizeof(long long) * _EventCodes.size() );
	}

    return retval;
}

int PapiProfiler::startThreadSet(papiThreadSet *thrset) {

    int retval;

    retval = PAPI_start(thrset->EventSet);
    if(retval != PAPI_OK) {
        handleError("Thread PAPI: Could not start event", retval);
    }

    return retval;
}

int PapiProfiler::stopThreadSet(papiThreadSet *thrset) {

    int retval;
    int *EventSet = &thrset->EventSet;
    long long *values = thrset->values; 
    retval = PAPI_stop( *EventSet, values);
    if ( retval != PAPI_OK ) {
        handleError("Thread PAPI: Could not stop events", retval);
    }

    return retval;
}

int PapiProfiler::clearThreadSet(papiThreadSet *thrset) {
    
    int retval;

    retval = PAPI_cleanup_eventset(thrset->EventSet);
    if(retval != PAPI_OK) {
        handleError("Thread PAPI: Could not clear event set", retval);
    }
    free(thrset->values);

    return retval;
}

int PapiProfiler::handleError(const char *msg, int retval) {
    printf("%s: %s\n", msg, PAPI_strerror(retval));
    return retval;
}

long long PapiProfiler::getEventCount(char *event) {

    long long values;
    return values;
}