#ifndef PAPIPROFILER_H_
#define PAPIPROFILER_H_

#ifdef __cplusplus
#include <vector>
#include <string>
#include <stdio.h> 
#include <stdlib.h>
#include "papi.h"
#include "log.h"
#endif

typedef struct _papiThreadSet{

	int EventSet;
	long long *values;

} papiThreadSet;

class PapiProfiler {

protected:

	std::vector<int> _EventCodes;
	papiThreadSet *_thrset_arr;
	int _num_threads;

public:

	PapiProfiler(int num_threads);
	~PapiProfiler();

	int init();
    int handleError(const char* msg, int retval);
    int addEvent(char *event);
    int clearEvents();
    long long getEventCount(char *event);

    void setNumThreads(int num_threads);
    int initThreadSet();
    int startThreadSet(papiThreadSet *thrset);
    int stopThreadSet(papiThreadSet *thrset);
    int clearThreadSet(papiThreadSet *thrset);

};


#endif /* PAPIPROFILER_H_ */