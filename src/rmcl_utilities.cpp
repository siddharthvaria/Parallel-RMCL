/*
 * rmcl_utilities.cpp
 *
 *  Created on: Mar 23, 2013
 *      Author: siddharth
 */

#include "rmcl.hpp"

void display_usage(){
	cout << "Please note:\n";
	cout << "Node indices start from zero\n";
	cout << "Input graph should be directed\n";
	cout << "[-f <number-input-files>] is no longer used\n";
	cout << "USAGE (Arguments can be in any order):\n";
	cout << "./rmcl_mpi [-e <number_edges(number of lines in input file)>] [-I <complete_input_filepath>] [-O <complete_output_filepath>] [-n <largest node index + 1>] [-w <n(unweighted) or y(weighted)>]\n";
	cout << "Optional:\n";
	cout << "[-c <converge threshold>]\t\tconverge threshold(default : 1e-5)\n";
	cout << "[-f <number-input-files>]\t\tnumber of input files to be read. If input graph is split into multiple files then -f is necessary(default : 1)\n";
	cout << "[-i <inflate-constant>]\t\tinflate constant(default : 2.0)\n";
	cout << "[-m <max-iterations>]\t\tmaximum number of iterations(default : 100)\n";
	cout << "[-h]\t\thelp\n";
	cout << "Example: ./rmcl_mpi -I /xyz/abc/pqr/input/filename -O /xyz/abc/pqr/output/filename -n 17000 -e 153000 -w n\n";
}

void printMatricesTriples(spMat *M,spMat *Mg,csi pid,csi nprocesses,csi nnp,csi nnp_exact){
	csi i;

	//print M
	cout << pid << ":Matrix M" << "(rows:" << M->nrows << ",columns:" << M->ncols << ")" <<  endl;
	cout << "(column,row,value):" << endl;
	for(i = 0;i < M->nz;i++){
		cout << "(" << M->cptrs[i] << "," << M->rinds[i] << "," << M->values[i] << ")" << " ";
	}

	cout << endl;

	//print Mg
	cout << pid << ":Matrix Mg" << "(rows:" << Mg->nrows << ",columns:" << Mg->ncols << ")" <<  endl;
	cout << "(column,row,value):" << endl;
	for(i = 0;i < Mg->nz;i++){
		cout << "(" << Mg->cptrs[i] << "," << Mg->rinds[i] << "," << Mg->values[i] << ")" << " ";
	}
	cout << endl;
}


void printMatricesCSC(spMat *M,spMat *Mg,csi pid,csi nprocesses,csi nnp,csi nnp_exact){
	csi i,j;

	//print M
	cout << pid << ":Matrix M" << "(rows:" << M->nrows << ",columns:" << M->ncols << ")" <<  endl;
	cout << "(column,row,value):" << endl;
	for(i = 0;i < M->ncols;i++){
		for(j = M->cptrs[i];j < M->cptrs[i+1];j++){
			cout << "(" << i << "," << M->rinds[j] << "," << M->values[j] << ")" << " ";
		}
	}

	cout << endl;

	//print Mg
	cout << pid << ":Matrix Mg" << "(rows:" << Mg->nrows << ",columns:" << Mg->ncols << ")" <<  endl;
	cout << "(column,row,value):" << endl;
	for(i = 0;i < Mg->ncols;i++){
		for(j = Mg->cptrs[i];j < Mg->cptrs[i+1];j++){
			cout << "(" << i << "," << Mg->rinds[j] << "," << Mg->values[j] << ")" << " ";
		}
	}
	cout << endl;
}

spMat** readDistribute(csi nnp,csi nnp_exact,csi nzmax,csi nzmax_exact,csi pid,csi nprocesses,long *dmc){

	ifstream inFile;
	csw colRowWeight[3];
	csi C,R;
	csi f,g;
	csw tmp;

	MPI_Status status;
	status.MPI_SOURCE = 0;
	status.MPI_TAG = 0;

	spMat** Matrices = (spMat**) calloc (2,sizeof(spMat*));
	if(Matrices == NULL){
		cout << "Matrices: Memory Allocation failed\n";
		MPI_Abort(MPI_COMM_WORLD,1);
	}
	*dmc += 2 * sizeof(spMat*);

	spMat *M,*Mg;

	if(pid == nprocesses - 1){
		M = spMat_spalloc (nnp_exact,arguments.nnodes,nzmax_exact,1,1,dmc);
		if(M == NULL){
			cout << "M: Memory Allocation failed\n";
			MPI_Abort(MPI_COMM_WORLD,1);
		}
		Mg = spMat_spalloc (arguments.nnodes,nnp_exact,nzmax_exact,1,1,dmc);
		if(Mg == NULL){
			cout << "Mg: Memory Allocation failed\n";
			MPI_Abort(MPI_COMM_WORLD,1);
		}
	}
	else{
		M = spMat_spalloc (nnp,arguments.nnodes,nzmax,1,1,dmc);
		if(M == NULL){
			cout << "M: Memory Allocation failed\n";
			MPI_Abort(MPI_COMM_WORLD,1);
		}
		Mg = spMat_spalloc (arguments.nnodes,nnp,nzmax,1,1,dmc);
		if(M == NULL){
			cout << "Mg: Memory Allocation failed\n";
			MPI_Abort(MPI_COMM_WORLD,1);
		}
	}

	if(pid == IOP){

		string line;
		colRowWeight[2] = 1;		//default weight : 1

			char filename[200];
			strcpy(filename,arguments.inputFileDirectory);
			inFile.open(filename);
			if(inFile.is_open()){
				while(getline(inFile,line) && line.length() != 0){
					if(line.at(0) == '#'){	//ignore comments
						continue;
					}
					stringstream crw(line);
					if(arguments.weightedOrNot == 'y' || arguments.weightedOrNot == 'Y'){
						crw >> colRowWeight[0] >> colRowWeight[1] >> colRowWeight[2];
					}
					else{
						crw >> colRowWeight[0] >> colRowWeight[1];
					}
					for(f = 0;f < 2;f++){
						C = (csi)colRowWeight[0] / nnp;
						R = (csi)colRowWeight[1] / nnp;

						if(C == IOP){
							if(!spMat_entry (Mg,(csi)colRowWeight[1],(csi)colRowWeight[0],colRowWeight[2],dmc)){
								cout << "line 160,spMat_entry Failed\n";
								MPI_Abort(MPI_COMM_WORLD,1);
							}
						}
						else{
							//mpi send
							MPI_Send(&colRowWeight,sizeof(colRowWeight),mpiw,C,0,MPI_COMM_WORLD);
						}
						if(R == IOP){
							if(!spMat_entry (M,(csi)colRowWeight[1],(csi)colRowWeight[0],colRowWeight[2],dmc)){
								cout << "line 170,spMat_entry Failed\n";
								MPI_Abort(MPI_COMM_WORLD,1);
							}
						}
						else if(R != C){
							//mpi send
							MPI_Send(&colRowWeight,sizeof(colRowWeight),mpiw,R,0,MPI_COMM_WORLD);
						}
						//swap colRowWeight[0] and colRowWeight[1]
						tmp = colRowWeight[0];
						colRowWeight[0] = colRowWeight[1];
						colRowWeight[1] = tmp;
					}//for
				}	//while
				inFile.close();
			}	//if
			else{
				cout << "Error opening input file!\n";
				MPI_Abort(MPI_COMM_WORLD,1);
			}

		colRowWeight[0] = colRowWeight[1] = colRowWeight[2] = 0;
		//change this to account for IOP
		for(f = 0;f < IOP;f++)
		MPI_Send(&colRowWeight,sizeof(colRowWeight),mpiw,f,1,MPI_COMM_WORLD);
		for(f = IOP + 1;f < nprocesses;f++)
		MPI_Send(&colRowWeight,sizeof(colRowWeight),mpiw,f,1,MPI_COMM_WORLD);
	}	//if
	else{
		while(1){
			MPI_Recv(&colRowWeight,sizeof(colRowWeight),mpiw,0,MPI_ANY_TAG,MPI_COMM_WORLD,&status);
			if(status.MPI_TAG)
				break;
			C = (csi)colRowWeight[0] / nnp;
			R = (csi)colRowWeight[1] / nnp;
			if(C == pid){
				if((colRowWeight[0] - (nnp * pid)) < 0 || !spMat_entry (Mg,(csi)colRowWeight[1],(csi)(colRowWeight[0] - (nnp * pid)),colRowWeight[2],dmc)){
					cout << "line 215,spMat_entry Failed\n";
					MPI_Abort(MPI_COMM_WORLD,1);
				}
			}
			if(R == pid){
				if((colRowWeight[1] - (nnp * pid)) < 0 || !spMat_entry (M,(csi)(colRowWeight[1] - (nnp * pid)),(csi)colRowWeight[0],colRowWeight[2],dmc)){
					cout << "line 226,spMat_entry Failed\n";
					MPI_Abort(MPI_COMM_WORLD,1);
				}
			}
		}	//while
	} //else

	//add self loops.
	if(pid == nprocesses - 1){
		for(f = 0;f < nnp_exact;f++){
			if(!spMat_entry (M,f,(nnp * pid) + f,1,dmc)){
				cout << "line 239,spMat_entry Failed\n";
				MPI_Abort(MPI_COMM_WORLD,1);
			}
			if(!spMat_entry (Mg,(nnp * pid) + f,f,1,dmc)){
				cout << "line 243,spMat_entry Failed\n";
				MPI_Abort(MPI_COMM_WORLD,1);
			}
		}
	}
	else{
		for(f = 0;f < nnp;f++){
			if(!spMat_entry (M,f,(nnp * pid) + f,1,dmc)){
				cout << "line 251,spMat_entry Failed\n";
				MPI_Abort(MPI_COMM_WORLD,1);
			}
			if(!spMat_entry (Mg,(nnp * pid) + f,f,1,dmc)){
				cout << "line 255,spMat_entry Failed\n";
				MPI_Abort(MPI_COMM_WORLD,1);
			}
		}
	}

	//convert M and Mg from triples to csc formats
	Matrices[1] = spMat_compress (M,dmc);
	Matrices[0] = spMat_compress (Mg,dmc);

	//free triples matrices
	spMat_spfree(M,dmc);
	spMat_spfree(Mg,dmc);

	//convert M and Mg into column stochastic
	csw *sumOfColumns = (csw *)calloc(SPMAT_MAX(arguments.nnodes,1),sizeof(csw));

	if(sumOfColumns == NULL){
		cout << "sumOfColumns: Memory Allocation failed\n";
		MPI_Abort(MPI_COMM_WORLD,1);
	}

	for(f = 0;f < Matrices[0]->ncols;f++){
		//one pass to get the sum
		for(g = Matrices[0]->cptrs[f];g < Matrices[0]->cptrs[f+1];g++){
			sumOfColumns[nnp * pid + f] += Matrices[0]->values[g];
		}
		//one pass to normalize the values
		for(g = Matrices[0]->cptrs[f];g < Matrices[0]->cptrs[f+1];g++){
			Matrices[0]->values[g] = Matrices[0]->values[g] / sumOfColumns[nnp * pid + f];
		}
	}

	// mpi all to all communication

	//csi* sendCounts = new (nothrow)csi[nprocesses];
	csi* sendCounts = (csi*)calloc(SPMAT_MAX(nprocesses,1),sizeof(csi));
	if(sendCounts == NULL){
		cout << "sentCounts: Memory Allocation failed\n";
		MPI_Abort(MPI_COMM_WORLD,1);
	}

	//csi* recvCounts = new (nothrow)csi[nprocesses];
	csi* recvCounts = (csi*) calloc(SPMAT_MAX(nprocesses,1),sizeof(csi));
	if(recvCounts == NULL){
		cout << "recvCounts: Memory Allocation failed\n";
		MPI_Abort(MPI_COMM_WORLD,1);
	}

	//csi* sendDisps = new (nothrow)csi[nprocesses];
	csi* sendDisps = (csi*) calloc(SPMAT_MAX(nprocesses,1),sizeof(csi));
	if(sendDisps == NULL){
		cout << "sentDisps: Memory Allocation failed\n";
		MPI_Abort(MPI_COMM_WORLD,1);
	}

	//csi* recvDisps = new (nothrow)csi[nprocesses];
	csi* recvDisps = (csi*) calloc(SPMAT_MAX(nprocesses,1),sizeof(csi));
	if(recvDisps == NULL){
		cout << "recvDisps: Memory Allocation failed\n";
		MPI_Abort(MPI_COMM_WORLD,1);
	}

	if(pid == nprocesses - 1){
		for(g = 0;g < nprocesses;g++){
			sendCounts[g] = nnp_exact;
		}
	}
	else{
		for(g = 0;g < nprocesses;g++){
			sendCounts[g] = nnp;
		}
	}

	sendCounts[pid] = 0;

	for(g = 0;g < nprocesses;g++){
		recvCounts[g] = nnp;
		sendDisps[g] = pid * nnp;
		recvDisps[g] = g * nnp;
	}

	recvCounts[nprocesses - 1] = nnp_exact;
	recvCounts[pid] = 0;

	MPI_Alltoallv(sumOfColumns,sendCounts,sendDisps,mpiw,sumOfColumns,recvCounts,recvDisps,mpiw,MPI_COMM_WORLD);

	MPI_Barrier(MPI_COMM_WORLD);

	// normalize Matrices[1]
	for(f = 0;f < arguments.nnodes;f++){
		for(g = Matrices[1]->cptrs[f];g < Matrices[1]->cptrs[f+1];g++){
			Matrices[1]->values[g] = Matrices[1]->values[g] / sumOfColumns[f];
		}
	}

	//free up allocated memory
	if(sumOfColumns) free(sumOfColumns);
	if(sendCounts) free(sendCounts);
	if(sendDisps) free(sendDisps);
	if(recvCounts) free(recvCounts);
	if(recvDisps) free(recvDisps);
	return Matrices;
}

void parallelMultiply(spMat **Matrices,csi nprocesses,csi pid,csi nnp,csi nnp_exact,csi *mgid,csi *nnz,csi maxnnz,long *dmc,csi iter){

	csi i;
	
	spMat **subMatrices = (spMat**) calloc(SPMAT_MAX(nprocesses,1),sizeof(spMat*));
	if(subMatrices == NULL){
		cout << "subMatrices: Memory Allocation failed\n";
		MPI_Abort(MPI_COMM_WORLD,1);
	}

	//get the first sub matrix before going inside the for loop
	if(mgid[pid] == nprocesses - 1){
		Matrices[0]->ncols = nnp_exact;
	}
	subMatrices[mgid[pid]] = spMat_multiply(Matrices[1],Matrices[0],dmc);
	Matrices[0]->ncols = nnp;

	spMat *newMg,*temp;
	csi currMgId;
	csi newMgId;
	csi sendPid;
	csi recvPid;
	csi sendSize;
	csi recvSize;
	MPI_Status status;
	
	//allocate newMg
	newMg = spMat_spalloc(arguments.nnodes,nnp,maxnnz,1,0,dmc);
	if(newMg == NULL){
		cout << "newMg: Memory Allocation failed\n";
		MPI_Abort(MPI_COMM_WORLD,1);
	}
	
	MPI_Barrier(MPI_COMM_WORLD);

	for(i = 1;i < nprocesses;i++){

		currMgId = mgid[(pid + i - 1) % nprocesses];
		newMgId = mgid[(pid+i) % nprocesses];
		sendPid = (pid > 0) ? pid - 1 : nprocesses - 1;
		recvPid = (pid + 1) % nprocesses;
		sendSize = (currMgId == (nprocesses - 1)) ? nnp_exact + 1 : nnp + 1;
		recvSize = (newMgId == (nprocesses - 1)) ? nnp_exact + 1 : nnp + 1;
/*		
		####################################
		changes made here
		newMg = spMat_spalloc(arguments.nnodes,recvSize - 1,nnz[newMgId],1,0,dmc);
		####################################
		
		if(newMg == NULL){
			cout << "newMg: Memory Allocation failed\n";
			MPI_Abort(MPI_COMM_WORLD,1);
		}
*/
		//send and recv cptrs

		MPI_Sendrecv(Matrices[0]->cptrs,sendSize,mpii,sendPid,\
		0,newMg->cptrs,recvSize,mpii,recvPid,0,MPI_COMM_WORLD,&status);

		MPI_Barrier(MPI_COMM_WORLD);

		//send and recv rinds
		MPI_Sendrecv(Matrices[0]->rinds,nnz[currMgId],mpii,sendPid,\
		0,newMg->rinds,nnz[newMgId],mpii,recvPid,0,MPI_COMM_WORLD,&status);

		MPI_Barrier(MPI_COMM_WORLD);

		//send and recv values
		MPI_Sendrecv(Matrices[0]->values,nnz[currMgId],mpiw,sendPid,\
		0,newMg->values,nnz[newMgId],mpiw,recvPid,0,MPI_COMM_WORLD,&status);

		MPI_Barrier(MPI_COMM_WORLD);
		
		//swap Matrices[0] and newMg
		temp = Matrices[0];
		Matrices[0] = newMg;
		newMg = temp;

		if(newMgId == nprocesses - 1){
			Matrices[0]->ncols = nnp_exact;
		}
		subMatrices[newMgId] = spMat_multiply(Matrices[1],Matrices[0],dmc);
		Matrices[0]->ncols = nnp;
	}

	//free newMg which is actually pointing to old block
	spMat_spfree(newMg,dmc);

	//merge submatrices and assign new M
	Matrices[1] = spMat_merge(subMatrices,nprocesses,dmc);
	if(subMatrices) free(subMatrices);
}

csw inflateAndPrune(spMat **Matrices,spMat *oldM,csi pid,csi nprocesses,long *dmc){
	int i,j;
	csw change;
	csw *partialSums = (csw *) calloc(SPMAT_MAX(arguments.nnodes,1),sizeof(csw));
	csw *globalChange;
	if(partialSums == NULL){
		cout << "partialSums: Memory Allocation failed\n";
		MPI_Abort(MPI_COMM_WORLD,1);
	}
/*
	old approach
	//inflate
	for(i = 0;i < arguments.nnodes;i++){
		for(j = Matrices[1]->cptrs[i];j < Matrices[1]->cptrs[i + 1];j++){
			//Matrices[1]->values[j] = pow((Matrices[1]->values[j] / partialSums[i]),arguments.inflateConstant);
			Matrices[1]->values[j] = pow((Matrices[1]->values[j]),arguments.inflateConstant);
			if(Matrices[1]->values[j] > 1e-6){
				partialSums[i] += Matrices[1]->values[j];
			}
			else{
				Matrices[1]->values[j] = 0;
			}
		}
	}

	//eliminate zero values from the matrix Matrices[1]
	spMat_dropzeros(Matrices[1]);

	//mpi all reduce in place
	MPI_Allreduce(MPI_IN_PLACE,partialSums,arguments.nnodes,mpiw,MPI_SUM,MPI_COMM_WORLD);

	for(i = 0;i < arguments.nnodes;i++){
		for(j = Matrices[1]->cptrs[i];j < Matrices[1]->cptrs[i + 1];j++){
			Matrices[1]->values[j] = Matrices[1]->values[j] / partialSums[i];
		}
	}

 */
	csi* nnz = (csi *) calloc(SPMAT_MAX(arguments.nnodes,1),sizeof(csi));
	if(nnz == NULL){
		cout << "nnz: Memory Allocation failed\n";
		MPI_Abort(MPI_COMM_WORLD,1);
	}

	csw* sum = (csw *) calloc(SPMAT_MAX(arguments.nnodes,1),sizeof(csw));
	if(sum == NULL){
		cout << "sum: Memory Allocation failed\n";
		MPI_Abort(MPI_COMM_WORLD,1);
	}

	csw* max = (csw *) malloc(SPMAT_MAX(arguments.nnodes,1) * sizeof(csw));
	if(max == NULL){
		cout << "max: Memory Allocation failed\n";
		MPI_Abort(MPI_COMM_WORLD,1);
	}

	csw* pruneThresholds;
	//inflate
	for(i = 0;i < arguments.nnodes;i++){
		max[i] = mincsw;
		for(j = Matrices[1]->cptrs[i];j < Matrices[1]->cptrs[i + 1];j++){
			//Matrices[1]->values[j] = pow((Matrices[1]->values[j] / partialSums[i]),arguments.inflateConstant);
			Matrices[1]->values[j] = pow((Matrices[1]->values[j]),arguments.inflateConstant);
			if(Matrices[1]->values[j] > 1e-6){
				sum[i] += Matrices[1]->values[j];
				nnz[i]++;
				if(Matrices[1]->values[j] > max[i]){
					max[i] = Matrices[1]->values[j];
				}
			}
			else{
				Matrices[1]->values[j] = 0;
			}
		}
	}

	pruneThresholds = computeThresholds(sum,nnz,max,pid);

	//free memories
	if(sum) free(sum);
	if(nnz) free(nnz);
	if(max) free(max);

	//prune and re-normalize
	for(i = 0;i < arguments.nnodes;i++){
		for(j = Matrices[1]->cptrs[i];j < Matrices[1]->cptrs[i + 1];j++){
			if(Matrices[1]->values[j] > pruneThresholds[i]){
				partialSums[i] += Matrices[1]->values[j];
			}
			else{
				Matrices[1]->values[j] = 0;
			}
		}
	}

	//free pruneThresholds
	if(pruneThresholds) free(pruneThresholds);

	//eliminate zero values from the matrix Matrices[1]
	spMat_dropzeros(Matrices[1],dmc);

	//mpi all reduce in place
	MPI_Allreduce(MPI_IN_PLACE,partialSums,arguments.nnodes,mpiw,MPI_SUM,MPI_COMM_WORLD);

	for(i = 0;i < arguments.nnodes;i++){
		for(j = Matrices[1]->cptrs[i];j < Matrices[1]->cptrs[i + 1];j++){
			Matrices[1]->values[j] = Matrices[1]->values[j] / partialSums[i];
		}
	}

	//compute change per column. reusing partialSums
	spMat_difference(Matrices[1],oldM,1,-1,partialSums);

	//free old Matrix M
	spMat_spfree(oldM,dmc);

	if(pid == IOP){
		globalChange = (csw*) malloc(SPMAT_MAX(arguments.nnodes,1) * sizeof(csw));
		if(globalChange == NULL){
			cout << "globalChange: Memory Allocation failed\n";
			MPI_Abort(MPI_COMM_WORLD,1);
		}
	}
	else{
		globalChange = NULL;
	}
	//gather global change at process IOP
	MPI_Reduce(partialSums,globalChange,arguments.nnodes,mpiw,MPI_SUM,IOP,MPI_COMM_WORLD);

	//free partialSums
	if(partialSums) free(partialSums);

	MPI_Barrier(MPI_COMM_WORLD);

	change = 0.0;
	if(pid == IOP){
		for(i = 0;i < arguments.nnodes;i++){
			change += sqrt(globalChange[i]);
		}

		//free globalChange
		if(globalChange) free(globalChange);

		change = change / arguments.nnodes;

		for(i = 1;i < nprocesses;i++){
			MPI_Send(&change,1,mpiw,i,0,MPI_COMM_WORLD);
		}
	}
	else{
		MPI_Status status;
		MPI_Recv(&change,1,mpiw,IOP,0,MPI_COMM_WORLD,&status);
	}

	MPI_Barrier(MPI_COMM_WORLD);

	return change;
}

void get_mgNonzeros(spMat *Mg,csi *nnz,csi *maxnnz,csi nprocesses,csi pid,csi nnp,csi nnp_exact){

	int i;

	nnz[pid] = Mg->cptrs[Mg->ncols];

	//mpi all to all communication

	//csi* sendCounts = new (nothrow)csi[nprocesses];
	csi* sendCounts = (csi*) calloc(SPMAT_MAX(nprocesses,1),sizeof(csi));
	if(sendCounts == NULL){
		cout << "sentCounts: Memory Allocation failed\n";
		MPI_Abort(MPI_COMM_WORLD,1);
	}

	//csi* recvCounts = new (nothrow)csi[nprocesses];
	csi* recvCounts = (csi*) calloc(SPMAT_MAX(nprocesses,1),sizeof(csi));
	if(recvCounts == NULL){
		cout << "recvCounts: Memory Allocation failed\n";
		MPI_Abort(MPI_COMM_WORLD,1);
	}

	//csi* sendDisps = new (nothrow)csi[nprocesses];
	csi* sendDisps = (csi*) calloc(SPMAT_MAX(nprocesses,1),sizeof(csi));
	if(sendDisps == NULL){
		cout << "sentDisps: Memory Allocation failed\n";
		MPI_Abort(MPI_COMM_WORLD,1);
	}

	//csi* recvDisps = new (nothrow)csi[nprocesses];
	csi* recvDisps = (csi*) calloc(SPMAT_MAX(nprocesses,1),sizeof(csi));
	if(recvDisps == NULL){
		cout << "recvDisps: Memory Allocation failed\n";
		MPI_Abort(MPI_COMM_WORLD,1);
	}


	for(i = 0;i < nprocesses;i++){
		sendCounts[i] = 1;
		recvCounts[i] = 1;
		sendDisps[i] = pid;
		recvDisps[i] = i;
	}

	sendCounts[pid] = 0;
	recvCounts[pid] = 0;

	MPI_Alltoallv(nnz,sendCounts,sendDisps,mpiw,nnz,recvCounts,recvDisps,mpiw,MPI_COMM_WORLD);

	//free allocated memory
	if(sendCounts) free(sendCounts);
	if(recvCounts) free(recvCounts);
	if(sendDisps) free(sendDisps);
	if(recvDisps) free(recvDisps);
	
	//each process computes the max size. One can compute and broadcast to others as well.
	
	for(i = 0;i < nprocesses;i++){
		if(nnz[i] > *maxnnz){
			*maxnnz = nnz[i];
		}
	}
}

csw* computeThresholds(csw* sum,csi* nnz,csw* max,csi pid){
	//
	csi i;
	csw* globalSum;
	csi* globalnnz;
	csw* globalMax;

	csw* thresholds = (csw*) calloc(SPMAT_MAX(arguments.nnodes,1),sizeof(csw));
	if(thresholds == NULL){
		cout << "thresholds: Memory Allocation failed\n";
		MPI_Abort(MPI_COMM_WORLD,1);
	}

	if(pid == IOP){
		globalSum = (csw*) malloc(SPMAT_MAX(arguments.nnodes,1) * sizeof(csw));
		if(globalSum == NULL){
			cout << "globalSum: Memory Allocation failed\n";
			MPI_Abort(MPI_COMM_WORLD,1);
		}

		globalnnz = (csi*) malloc(SPMAT_MAX(arguments.nnodes,1) * sizeof(csi));
		if(globalnnz == NULL){
			cout << "globalnnz: Memory Allocation failed\n";
			MPI_Abort(MPI_COMM_WORLD,1);
		}

		globalMax = (csw*) malloc(SPMAT_MAX(arguments.nnodes,1) * sizeof(csw));
		if(globalMax == NULL){
			cout << "globalMax: Memory Allocation failed\n";
			MPI_Abort(MPI_COMM_WORLD,1);
		}
	}
	else{
		globalSum = NULL;
		globalnnz = NULL;
		globalMax = NULL;
	}

	//gather globalSum at IOP
	MPI_Reduce(sum,globalSum,arguments.nnodes,mpiw,MPI_SUM,IOP,MPI_COMM_WORLD);
	//gather globalnnz at IOP
	MPI_Reduce(nnz,globalnnz,arguments.nnodes,mpii,MPI_SUM,IOP,MPI_COMM_WORLD);
	//gather max at IOP
	MPI_Reduce(max,globalMax,arguments.nnodes,mpiw,MPI_MAX,IOP,MPI_COMM_WORLD);

	MPI_Barrier(MPI_COMM_WORLD);

	if(pid == IOP){
		csw avg;
		for(i = 0;i < arguments.nnodes;i++){
			avg = globalSum[i] / globalnnz[i];
			thresholds[i] = RMCL_PRUNE_A * avg * (1 - RMCL_PRUNE_B * (globalMax[i]-avg));
		    thresholds[i] = (thresholds[i] > globalMax[i]) ? globalMax[i] : abs(thresholds[i]);
		}

		//free memories
		if(globalSum) free(globalSum);
		if(globalnnz) free(globalnnz);
		if(globalMax) free(globalMax);
	}

	//mpi all reduce in place. All but IOP's thresholds are zero
	MPI_Allreduce(MPI_IN_PLACE,thresholds,arguments.nnodes,mpiw,MPI_SUM,MPI_COMM_WORLD);

	return thresholds;
}

void interpretClusters(spMat** Matrices,csi pid,csi nnp){

	pairs* partialMaxAndIndices =  (pairs*) malloc(SPMAT_MAX(arguments.nnodes,1) * sizeof(pairs));
	if(partialMaxAndIndices == NULL){
		cout << "partialMaxAndIndices: Memory Allocation failed\n";
		MPI_Abort(MPI_COMM_WORLD,1);
	}
	pairs* globalMaxAndIndices;
	csi i,j;
	for(i = 0;i < arguments.nnodes;i++){
		partialMaxAndIndices[i].mxval = mincsw;
		for(j = Matrices[1]->cptrs[i];j < Matrices[1]->cptrs[i+1];j++){
			if(Matrices[1]->values[j] > partialMaxAndIndices[i].mxval){
				partialMaxAndIndices[i].mxval = Matrices[1]->values[j];
				partialMaxAndIndices[i].ind = pid * nnp + Matrices[1]->rinds[j];
			}
		}
	}

	if(pid == IOP){
		globalMaxAndIndices = (pairs*) malloc(SPMAT_MAX(arguments.nnodes,1) * sizeof(pairs));
		if(globalMaxAndIndices == NULL){
			cout << "globalMaxAndIndices: Memory Allocation failed\n";
			MPI_Abort(MPI_COMM_WORLD,1);
		}
	}
	else{
		globalMaxAndIndices = NULL;
	}

	MPI_Reduce(partialMaxAndIndices,globalMaxAndIndices,arguments.nnodes,mpiwi,MPI_MAXLOC,IOP,MPI_COMM_WORLD);

	if(partialMaxAndIndices) free(partialMaxAndIndices);

	if(pid == IOP){
		ofstream outFile;

		outFile.open(arguments.outputFileDirectory);

		if(outFile.is_open()){
			for(i = 0;i < arguments.nnodes;i++){
				outFile << globalMaxAndIndices[i].ind << endl;
			}
			outFile.close();
		}
		else{
			cout << "Error creating output file!\n";
			MPI_Abort(MPI_COMM_WORLD,1);
		}
		if(globalMaxAndIndices) free(globalMaxAndIndices);
	}
}

void resizeMgBlocks(spMat *Mg,csi maxnnz,csi pid,csi nnp,csi nprocesses,long *dmc){
	csi okj;
	if(pid == nprocesses - 1){
		//resize cptrs
		Mg->cptrs = (csi *) spMat_realloc (Mg->cptrs,nnp + 1, sizeof (csi), &okj);
		if(!okj){
			cout << "MgBlock: Memory Reallocation failed\n";
			MPI_Abort(MPI_COMM_WORLD,1);
		}
		*dmc += (nnp - Mg->ncols) * sizeof(csi);
		Mg->ncols = nnp;
	}

	//resize rinds and values
	if(Mg->nzmax < maxnnz){
		if(!spMat_sprealloc (Mg,maxnnz,dmc)){
			cout << "Matrices[0]: Memory Reallocation failed\n";
		}
	}
}


