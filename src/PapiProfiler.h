#ifndef PAPIPROFILER_H_
#define PAPIPROFILER_H_

#ifdef __cplusplus
#include <vector>
#include <string>
#include <stdio.h> 
#include <stdlib.h>
#include "papi.h"
#include "log.h"

#define CHARBUF 128
// #define THR_SEPARATE -1
// #define THR_REDUCE -2
// #define THR_SEPARATE_GBL -3
// #define THR_REDUCE_GBL -4

#endif

enum {
	THR_SEPARATE,
	THR_REDUCE,
	THR_SEPARATE_GBL,
	THR_REDUCE_GBL
};

typedef struct _papiThreadSet{

	int EventSet;
	long long *values;
	bool init;

} papiThreadSet;

class PapiProfiler {

protected:

	std::vector< std::vector<int> > _CodeSections;

	/* Vector of the codes of Events to count */
	std::vector<int> _EventCodes;

	/* Vector of papiThreadSet structs which keep track of
	   counts for each thread, of length _NUM_THREADS */	
	std::vector<papiThreadSet> _thrset_arr;
	std::vector<papiThreadSet> _gbl_thrset_arr;

	int _num_codeSections;
	int _num_threads;

public:

	PapiProfiler(int num_threads, int num_codeSections);
	~PapiProfiler();

	int init();
    int handleError(const char* msg, int retval);
    int addEvent(char *event);
    int getNumEvents();
    int clearEvents();
    void printEventCountsPerUnit(int reduce, int perunit);
    void printEventCounts(int reduce);
    void gblAccum();

    void setNumThreads(int num_threads);
    int initThreadSet(int tid);
    int startThreadSet(int tid);
    int accumThreadSet(int tid);
    int stopThreadSet(int tid);
    int clearThreadSet();

};


#endif /* PAPIPROFILER_H_ */