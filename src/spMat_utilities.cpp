/*
 * spMat_utilities.cpp
 *
 *  Created on: Mar 23, 2013
 *      Author: siddharth
 */

/*
CSparse: a Concise Sparse matrix package.
Copyright (c) 2006, Timothy A. Davis.
http://www.cise.ufl.edu/research/sparse/CSparse

--------------------------------------------------------------------------------

CSparse is free software; you can redistribute it and/or
modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation; either
version 2.1 of the License, or (at your option) any later version.

CSparse is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
Lesser General Public License for more details.
*/

#include "rmcl.hpp"
//#include "../Include/spMatrix.h"

/*
cs = spMat
m = nrows
n = ncols
p = cptrs
i = rinds
x = values
*/

/*
void* spMat_calloc (csi n, size_t size) {
    return (calloc (SPMAT_MAX (n,1), size)) ;
}

void *spMat_malloc (csi n, size_t size)
{
    return (malloc (SPMAT_MAX (n,1) * size)) ;
}

void *spMat_free (void *p)
{
    if (p) free (p) ;       // free p if it is not already NULL
    return (NULL) ;         //return NULL to simplify the use of spMat_free
}
*/
 void *spMat_realloc (void *p, csi n, size_t size, csi *ok)
{
    void *pnew ;
    pnew = realloc (p, SPMAT_MAX (n,1) * size) ; /* realloc the block */
    *ok = (pnew != NULL) ;                  /* realloc fails if pnew is NULL */
    return ((*ok) ? pnew : p) ;             /* return original p if failure */
}

spMat *spMat_spalloc (csi nrows, csi ncols, csi nzmax, csi values, csi triplet,long *dmc)
{
    //spMat *A = (spMat *)spMat_calloc (1, sizeof(spMat)) ;    /* allocate the spMat struct */
	spMat *A = (spMat *)calloc (1,sizeof(spMat));
    if (!A) return (NULL) ;                 /* out of memory */
	*dmc += 1 * sizeof(spMat);
    A->nrows = nrows ;                              /* define dimensions and nzmax */
    A->ncols = ncols ;
    A->nzmax = nzmax = SPMAT_MAX (nzmax, 1) ;
    A->nz = triplet ? 0 : -1 ;              /* allocate triplet or comp.col. triplet:0 for compressed sparse columns */
    //A->cptrs = (csi*) spMat_malloc (triplet ? nzmax : ncols+1, sizeof (csi)) ;
    A->cptrs = (csi*) malloc(SPMAT_MAX(triplet ? nzmax : ncols+1,1) * sizeof (csi));
    *dmc += SPMAT_MAX(triplet ? nzmax : ncols+1,1) * sizeof (csi);
    //A->rinds = (csi*) spMat_malloc (nzmax, sizeof (csi)) ;
    A->rinds = (csi*) malloc(SPMAT_MAX(nzmax,1) * sizeof (csi));
    *dmc += SPMAT_MAX(nzmax,1) * sizeof (csi);
    //A->values = values ? (csw *) spMat_malloc (nzmax,sizeof (csw)) : NULL ;
    A->values = values ? (csw *) malloc (SPMAT_MAX(nzmax,1) * sizeof (csw)) : NULL ;
    *dmc += values ? SPMAT_MAX(nzmax,1) * sizeof (csw) : 0;
    return ((!A->cptrs || !A->rinds || (values && !A->values)) ? spMat_spfree (A,dmc) : A) ;
}

spMat *spMat_spfree (spMat *A,long *dmc)
{
    if (!A) return (NULL) ;     /* do nothing if A already NULL */
    //spMat_free (A->cptrs) ;
    if (A->cptrs) {
    	*dmc -= SPMAT_MAX(SPMAT_CSC(A) ? A->ncols+1 : A->nzmax,1) * sizeof (csi);
    	free (A->cptrs);
    }
    //spMat_free (A->rinds) ;
    if (A->rinds) {
    	*dmc -= SPMAT_MAX(A->nzmax,1) * sizeof (csi);
    	free (A->rinds);
    }
    //spMat_free (A->values) ;
    if (A->values) {
    	*dmc -= SPMAT_MAX(A->nzmax,1) * sizeof (csw);
    	free (A->values);
    }
    //return ((spMat*) spMat_free (A)) ;      /* free the spMat struct and return NULL */
    if (A) {
    	*dmc -= 1 * sizeof(spMat);
    	free (A) ;       /* free p if it is not already NULL */
    }
    return ((spMat*)NULL) ;
}

//increases the maximum number of entries the matrix can hold
csi spMat_sprealloc (spMat *A, csi nzmax,long *dmc)
{
    csi ok, okrinds, okj = 1, okvalues = 1 ;
    if (!A) return (0) ;
    if (nzmax <= 0) nzmax = (SPMAT_CSC (A)) ? (A->cptrs [A->ncols]) : A->nz ;
    A->rinds = (csi*) spMat_realloc (A->rinds, nzmax, sizeof (csi), &okrinds) ;
    if(okrinds){
    	*dmc += ((SPMAT_MAX (nzmax,1) - A->nzmax) * sizeof(csi));
    }
    if (SPMAT_TRIPLET (A)) {
    	A->cptrs = (csi*) spMat_realloc (A->cptrs, nzmax, sizeof (csi), &okj) ;
    	if(okj){
    		*dmc += ((SPMAT_MAX (nzmax,1) - A->nzmax) * sizeof(csi));
    	}
    }
    if (A->values) {
    	A->values = (csw *) spMat_realloc (A->values, nzmax, sizeof (csw), &okvalues) ;
    	if(okvalues){
    		*dmc += ((SPMAT_MAX (nzmax,1) - A->nzmax) * sizeof(csw));
    	}
    }
    ok = (okrinds && okj && okvalues);
    if (ok) {
    	A->nzmax = nzmax;
    }
    else{
    	if(!okrinds)
        	*dmc -= ((SPMAT_MAX (nzmax,1) - A->nzmax) * sizeof(csi));
    	if(!okj)
        	*dmc -= ((SPMAT_MAX (nzmax,1) - A->nzmax) * sizeof(csi));
    	if(!okvalues)
        	*dmc -= ((SPMAT_MAX (nzmax,1) - A->nzmax) * sizeof(csw));
    }
    return (ok) ;
}

/* add an entry to a triplet matrix; return 1 if ok, 0 otherwise */
csi spMat_entry (spMat *T, csi r, csi c, csw v,long *dmc)
{
    if (!T || (T->nz >= T->nzmax && !spMat_sprealloc (T, 2 * (T->nzmax),dmc))) return(0);
    if (T->values) T->values [T->nz] = v ;
    T->rinds [T->nz] = r ;
    T->cptrs [T->nz++] = c ;
    T->nrows = SPMAT_MAX (T->nrows, r+1) ;
    T->ncols = SPMAT_MAX (T->ncols, c+1) ;
    return (1) ;
}

spMat *spMat_compress (const spMat *T,long *dmc)
{
    csi m, n, nz, p, k, *Cp, *Ci, *w, *Ti, *Tj ;
    csw *Cx, *Tx ;
    spMat *C ;
    if (!SPMAT_TRIPLET (T)) return (NULL) ;                /* check inputs */
    m = T->nrows ; n = T->ncols ; Ti = T->rinds ; Tj = T->cptrs ; Tx = T->values ; nz = T->nz;
    C = spMat_spalloc (m, n, nz, Tx != NULL, 0,dmc) ;          /* allocate result. 0 for CSC format */
    //w = (csi *) spMat_calloc (n, sizeof (csi)) ;                   /* get workspace */
    w = (csi *) calloc (SPMAT_MAX(n,1), sizeof (csi)) ;
    if (!C || !w){
    	return (spMat_done (C, w, NULL,0,dmc)) ;    /* out of memory. 0 for failure */
    }
    Cp = C->cptrs ; Ci = C->rinds ; Cx = C->values ;
    for (k = 0 ; k < nz ; k++) w [Tj [k]]++ ;           /* column counts */
    spMat_cumsum (Cp, w, n) ;                              /* column pointers */
    for (k = 0 ; k < nz ; k++)
    {
        Ci [p = w [Tj [k]]++] = Ti [k] ;    /* A(i,j) is the pth entry in C */
        if (Cx) Cx [p] = Tx [k] ;
    }
    return (spMat_done (C, w, NULL, 1,dmc)) ;      /* success; free w and return C */
}

spMat *spMat_done (spMat *C, void *w, void *x, csi ok,long *dmc)
{
    //spMat_free (w) ;                       /* free workspace */
    if(w) free(w);
    //spMat_free (x) ;
    if(x) free(x);
    return (ok ? C : spMat_spfree (C,dmc)) ;   /* return result if OK, else free it */
}

double spMat_cumsum (csi *p, csi *c, csi n)
{
    csi i, nz = 0 ;
    double nz2 = 0 ;
    if (!p || !c) return (-1) ;     /* check inputs */
    for (i = 0 ; i < n ; i++)
    {
        p [i] = nz ;
        nz += c [i] ;
        nz2 += c [i] ;              /* also in double to avoid csi overflow */
        c [i] = p [i] ;             /* also copy p[0..n-1] back into c[0..n-1]*/
    }
    p [n] = nz ;
    return (nz2) ;                  /* return sum (c [0..n-1]) */
}
/*
csi spMat_nonzero (csi i, csi j, csw aij, void *other)
{
    return (aij != 0) ;
}
*/
csi spMat_dropzeros (spMat *A,long *dmc)
{
    //return (spMat_fkeep (A, &spMat_nonzero, NULL)) ;  /* keep all nonzero entries */
	return (spMat_fkeep (A,dmc)) ;  /* keep all nonzero entries */
}

csi spMat_fkeep (spMat *A,long *dmc)
{
    csi j, p, nz = 0, n, *Ap, *Ai ;
    csw *Ax ;
    if (!SPMAT_CSC (A)) return (-1) ;    /* check inputs */
    n = A->ncols ; Ap = A->cptrs ; Ai = A->rinds ; Ax = A->values ;
    for (j = 0 ; j < n ; j++)
    {
        p = Ap [j] ;                        /* get current location of col j */
        Ap [j] = nz ;                       /* record new location of col j */
        for ( ; p < Ap [j+1] ; p++)
        {
        	/*
            if (fkeep (Ai [p], j, Ax ? Ax [p] : 1, other))
            {
                if (Ax) Ax [nz] = Ax [p] ;  //keep A(i,j)
                Ai [nz++] = Ai [p] ;
            }
        	*/
            if ((Ax ? Ax [p] : 1) != 0)
            {
                if (Ax) Ax [nz] = Ax [p] ;  /* keep A(i,j) */
                Ai [nz++] = Ai [p] ;
            }
        }
    }
    Ap [n] = nz ;                           /* finalize A */
    if(!spMat_sprealloc (A, 0,dmc)) /* remove extra space from A */{
    	cout << "From spMat_fkeep: spMat_sprealloc failed\n";
    }
    return (nz) ;
}

spMat* spMat_merge(spMat** subMatrices,csi nprocesses,long *dmc){

	csi i,j;
	csi total_nnz = 0;
	csi disp = 0;
	csi partial_nnz = 0;

	for(i = 0;i < nprocesses;i++){
		total_nnz += subMatrices[i]->cptrs[subMatrices[i]->ncols];
	}

	spMat* newM = spMat_spalloc(subMatrices[0]->nrows,arguments.nnodes,total_nnz,1,0,dmc);
	if(newM == NULL){
		cout << "newM: Memory Allocation failed\n";
		MPI_Abort(MPI_COMM_WORLD,1);
	}

	newM->cptrs[0] = 0;

	for(i = 0;i < nprocesses;i++){
		partial_nnz = newM->cptrs[disp];
		for(j = 0;j <= subMatrices[i]->ncols;j++){
			newM->cptrs[disp + j] = subMatrices[i]->cptrs[j] + partial_nnz;
		}
		for(j = 0;j < subMatrices[i]->cptrs[subMatrices[i]->ncols];j++){
			newM->rinds[partial_nnz + j] = subMatrices[i]->rinds[j];
			newM->values[partial_nnz + j] = subMatrices[i]->values[j];
		}
		disp += subMatrices[i]->ncols;
		spMat_spfree(subMatrices[i],dmc);
	}

	return newM;
}

void spMat_difference (const spMat *A, const spMat *B, csw alpha, csw beta,csw * change)
{
    csi p, j, nz;
    //csi anz;
    //csi *Bp;
    csi m, n;
    //csi bnz;
    csi *w, values;
    csw *x, *Bx;
    //spMat *C  = (spMat*) spMat_malloc(1,sizeof(spMat));
    spMat *C  = (spMat*) calloc(1,sizeof(spMat));
    C->nz = -1;		//for CSC
    if (!SPMAT_CSC (A) || !SPMAT_CSC (B)) return;         /* check inputs */
    if (A->nrows != B->nrows || A->ncols != B->ncols) return;
    m = A->nrows ;
    //anz = A->cptrs [A->ncols] ;
    n = B->ncols ;
    //Bp = B->cptrs ;
    Bx = B->values ;
    w = (csi*) calloc(SPMAT_MAX(m,1), sizeof (csi)) ;                       /* get workspace */
    values = (A->values != NULL) && (Bx != NULL) ;
    //x = values ? (csw*)spMat_malloc (m, sizeof (csw)) : NULL ;    /* get workspace */
    x = values ? (csw*) malloc(SPMAT_MAX(m,1) * sizeof (csw)) : NULL ;    /* get workspace */
    //C->rinds = (csi*) spMat_malloc (m, sizeof (csi));
    C->rinds = (csi*) malloc(SPMAT_MAX(m,1) * sizeof (csi));
    if (!C || !C->rinds || !w || (values && !x)){
    	cout << "spMat_difference:  Memory Allocation failed\n";
    	if(C->rinds) free(C->rinds);
    	if(C) free(C);
    	if(w) free(w);
    	if(x) free(x);
        return;
    }

    for (j = 0 ; j < n ; j++)
    {
    	change[j] = 0;
    	nz = 0;
        nz = spMat_scatter (A, j, alpha, w, x, j+1, C, nz) ;   /* alpha*A(:,j)*/
        nz = spMat_scatter (B, j, beta, w, x, j+1, C, nz) ;    /* beta*B(:,j) */
        if (values) for (p = 0; p < nz ; p++) change[j] += pow(x[C->rinds[p]],2);
    }

    if(C->rinds) free(C->rinds);
    if(C) free(C);
    if(w) free(w);
    if(x) free(x);
    return;
}
