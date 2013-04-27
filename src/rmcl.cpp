/*
 * rmcl.cpp
 *
 *  Created on: Mar 23, 2013
 *      Author: siddharth
 */

#include "rmcl.hpp"

commandLineArguments arguments;

	int main(int argc,char *argv[]){

		long *dmc;		//keeps track of allocated dynamic memory
		dmc = (long *)calloc(1,sizeof(long));
		csi opt = 0;
		const char *optStr = "c:e:f:i:I:m:n:O:w:h";

		csi pid;
		csi nprocesses;

		MPI_Init(&argc,&argv);
		MPI_Comm_rank(MPI_COMM_WORLD,&pid);
		MPI_Comm_size(MPI_COMM_WORLD,&nprocesses);

		//cout << argv[0];
		if(argc == 1 || (argc > 1 && (strcmp(argv[1],"-h") == 0))){
			if(pid == IOP){
				display_usage();
			}
			MPI_Finalize();
			exit(EXIT_SUCCESS);
		}
		else if(argc < 11){
			if(pid == IOP){
				cout << "Not enough arguments!\n";
				display_usage();
			}
			MPI_Finalize();
			exit(EXIT_FAILURE);
		}
		// Initialize arguments to default parameters.
		arguments.convergeThreshold = 1e-5;
		arguments.fileCount = 1;
		arguments.inflateConstant = 2.0;
		arguments.inputFileDirectory = NULL;
		arguments.maxIterations = 100;
		arguments.nnodes = 0;
		arguments.nedges = 0;
		arguments.outputFileDirectory = NULL;
		arguments.weightedOrNot = 'n';

		// Process the arguments with getopt().

		opt = getopt(argc,argv,optStr);

		while( opt != -1 ) {
			switch( opt ) {
				case 'c':
					arguments.convergeThreshold = atof(optarg);
					break;

				case 'e':
					arguments.nedges = atoi(optarg) * 2;
					break;

				case 'f':
					arguments.fileCount = atoi(optarg);
					break;

				case 'i':
					arguments.inflateConstant = atof(optarg);
					break;

				case 'I':
					arguments.inputFileDirectory = optarg;
					break;

				case 'm':
					arguments.maxIterations = atoi(optarg);
					break;

				case 'n':
					arguments.nnodes = atoi(optarg);
					break;

				case 'O':
					arguments.outputFileDirectory = optarg;
					break;

				case 'w':
					arguments.weightedOrNot = optarg[0];
					break;

				case 'h':	/* fall-through is intentional */
					if(pid == IOP){
					display_usage();
					}
					MPI_Finalize();
					exit(EXIT_SUCCESS);
					//MPI_Abort(MPI_COMM_WORLD,1);

				default:
					cout << "Unknown Argument found!\n";
					if(pid == IOP){
					display_usage();
					}
					MPI_Finalize();
					exit(EXIT_FAILURE);
					//MPI_Abort(MPI_COMM_WORLD,1);
			}

			opt = getopt( argc, argv, optStr );
		}

		ofstream outFile;

		csi nzmax,nzmax_exact;		// nzmax_exact used only by the last process
		csi nnp,nnp_exact;			// nnp_exact used only by the last process

		nnp = ((arguments.nnodes - 1) / nprocesses) + 1;
		//if there is a need to think about the number of nodes in the last process
		nnp_exact = arguments.nnodes - (nnp * (nprocesses - 1));

		//to account for self loops
		arguments.nedges += arguments.nnodes;

		//estimate initial nzmax to minimize memory reallocation
		nzmax = (ceil(((csw)arguments.nedges) / arguments.nnodes) * nnp);
		//if there is a need to think about nzmax in the last process
		nzmax_exact = (ceil(((csw)arguments.nedges) / arguments.nnodes) * nnp_exact);

#ifdef DEBUG
		if(pid == IOP){
			cout << "nnp:" << nnp << endl;
			cout << "nnp_exact:" << nnp_exact << endl;
			cout << "nzmax:" << nzmax << endl;
			cout << "nzmax_exact:" << nzmax_exact << endl;
		}
#endif

		MPI_Barrier(MPI_COMM_WORLD);

		//File I/O done by process IOP
		spMat **Matrices;

		Matrices = readDistribute(nnp,nnp_exact,nzmax,nzmax_exact,pid,nprocesses,dmc);

		spMat *M;
		M = Matrices[1];

		MPI_Barrier(MPI_COMM_WORLD);

		//start iterative computations

		csw change = 100;
		csi i,iter = 0;
		csi lstmgid;
		double totalTimeElapsed;
		double timeElapsed;

		//mgid keeps track of the ownership of Mg blocks
		//csi* mgid = new (nothrow)csi[nprocesses];
		csi* mgid = (csi*)calloc(SPMAT_MAX(nprocesses,1),sizeof(csi));
		if(mgid == NULL){
			cout << "mgid: Memory Allocation failed\n";
			MPI_Abort(MPI_COMM_WORLD,1);
		}

		for(i = 0;i < nprocesses;i++){
			mgid[i] = i;
		}

		//csi* mgnnz = new (nothrow)csi[nprocesses];
		csi* mgnnz = (csi*)calloc(SPMAT_MAX(nprocesses,1),sizeof(csi));
		if(mgnnz == NULL){
			cout << "mgnnz: Memory Allocation failed\n";
			MPI_Abort(MPI_COMM_WORLD,1);
		}


		csi max_mgnnz = 0;	//stores the max(Mg) size 
		get_mgNonzeros(Matrices[0],mgnnz,&max_mgnnz,nprocesses,pid,nnp,nnp_exact);
		resizeMgBlocks(Matrices[0],max_mgnnz,pid,nnp,nprocesses,dmc);

		MPI_Barrier(MPI_COMM_WORLD);
		
#ifdef DEBUG
		cout << "Starting Iterative computation\n";
#endif
		totalTimeElapsed = MPI_Wtime();
		while(change > arguments.convergeThreshold && iter < arguments.maxIterations){
			timeElapsed = MPI_Wtime();
			parallelMultiply(Matrices,nprocesses,pid,nnp,nnp_exact,mgid,mgnnz,max_mgnnz,dmc,iter);
			change = inflateAndPrune(Matrices,M,pid,nprocesses,dmc);
			//assign new M to M (as old M) for next iteration
			M = Matrices[1];

			lstmgid = mgid[nprocesses - 1];
			for(i = nprocesses - 1;i > 0;i--){
				mgid[i] = mgid[i - 1];
			}
			mgid[0] = lstmgid;
			iter++;
			if(pid == IOP){
				cout << "##################\t" << "End of iteration:" << iter << "\t##################\n";
				cout << "Time elapsed this iteration:" << (MPI_Wtime() - timeElapsed) << " seconds\n";
				cout << "Change:" << change << endl;
			}

		}

		interpretClusters(Matrices,pid,nnp);

		if(pid == IOP){
			cout << "##################   Overall   ##################\n";
			cout << "Total number of iterations:" << iter << endl;
			cout << "Total running time:" << (MPI_Wtime() - totalTimeElapsed) << " seconds\n";
		}

		//free matrices and all memory
		spMat_spfree(Matrices[0],dmc);
		spMat_spfree(Matrices[1],dmc);
		if(Matrices) free(Matrices);
		if(mgid) free(mgid);
		if(mgnnz) free(mgnnz);
		cout << pid << ":Dynamic memory allocated at termination:" << *dmc << " bytes\n";
		if(dmc) free(dmc);
		MPI_Finalize();
}

