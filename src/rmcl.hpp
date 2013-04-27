/*
 * rmcl.hpp
 *
 *  Created on: Mar 23, 2013
 *      Author: siddharth
 */

#ifndef RMCL_HPP_
#define RMCL_HPP_


#include <cstdlib>
#include <climits>
#include <cfloat>
#include <cmath>
#include <cstdio>
#include <iostream>
//#include <string>
#include <cstring>
#include <fstream>
#include <sstream>
#include <cstddef>
#include <unistd.h>
#include <mpi.h>
#include "spMat.hpp"

using namespace std;

//#define DEBUG			//uncomment to print some messages to stdout
#define IOP	0
#define RMCL_PRUNE_A	0.90	/* pruning parameter */
#define RMCL_PRUNE_B	2		/* pruning parameter */

typedef struct{
	csw mxval;
	csi ind;
} pairs;

struct commandLineArguments{
	float convergeThreshold;	/* used to check of rmcl convergence (-c) */
	csi fileCount;			/* number of input files to be read (-f)*/
	float inflateConstant;		/* constant used by rmcl (-i) */
	char * inputFileDirectory;	/* path upto the directory where the input file is present (-I) */
	csi maxIterations;			/* max iterations of rmcl (-m)*/
	csi nnodes;			/* number of nodes (-n) */
	char * outputFileDirectory; /* path upto the directory where the output file is expected (-O) */
	char weightedOrNot;		/* indicates whether the graph is weighted or unweighted */
	csi nedges;				/* number of edges in the graph */
};

extern commandLineArguments arguments;

void display_usage();

spMat** readDistribute(csi nnp,csi nnp_exact,csi nzmax,csi nzmax_exact,csi pid,csi nprocesses,long *dmc);

void printMatricesCSC(spMat *M,spMat *Mg,csi pid,csi nprocesses,csi nnp,csi nnp_exact);

void printMatricesTriples(spMat *M,spMat *Mg,csi pid,csi nprocesses,csi nnp,csi nnp_exact);

void parallelMultiply(spMat **Matrices,csi nprocesses,csi pid,csi nnp,csi nnp_exact,csi *mgid,csi *nnz,csi maxnnz,long *dmc,csi iter);

void get_mgNonzeros(spMat *Mg,csi * nnz,csi *maxnnz,csi nprocesses,csi pid,csi nnp,csi nnp_exact);

csw inflateAndPrune(spMat **Matrices,spMat *oldM,csi pid,csi nprocesses,long *dmc);

csw* computeThresholds(csw* sum,csi* nnz,csw* max,csi pid);

void interpretClusters(spMat** Matrices,csi pid,csi nnp);

void resizeMgBlocks(spMat *Mg,csi maxnnz,csi pid,csi nnp,csi nprocesses,long *dmc);

#endif /* RMCL_HPP_ */
