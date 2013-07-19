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
	bool init;

} papiThreadSet;

class PapiProfiler {

protected:

	std::vector<int> _EventCodes;
	std::vector<papiThreadSet> _thrset_arr;
	int _num_threads;

public:

	PapiProfiler(int num_threads);
	~PapiProfiler();

	int init();
    int handleError(const char* msg, int retval);
    int addEvent(char *event);
    int getNumEvents();
    int clearEvents();
    void printEventCounts();

    void setNumThreads(int num_threads);
    int initThreadSet(int tid);
    int startThreadSet(int tid);
    int accumThreadSet(int tid);
    int stopThreadSet(int tid);
    int clearThreadSet();

};


#endif /* PAPIPROFILER_H_ */