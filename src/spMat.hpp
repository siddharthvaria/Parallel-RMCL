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


/*
 * spMat.hpp
 *
 *  Created on: Mar 23, 2013
 *      Author: siddharth
 */


#ifndef SPMAT_HPP_
#define SPMAT_HPP_


#define csi int
#define mpii MPI_INT

#define csw float
#define mpiw MPI_FLOAT
#define mpiwi MPI_FLOAT_INT
#define mincsw FLT_MIN



#define SPMAT_MAX(a,b) (((a) > (b)) ? (a) : (b))
#define SPMAT_MIN(a,b) (((a) < (b)) ? (a) : (b))
#define SPMAT_CSC(A) (A && (A->nz == -1))
#define SPMAT_TRIPLET(A) (A && (A->nz >= 0))

/*
cs = spMat
m = nrows
n = ncols
p = cptrs
i = rinds
x = values
*/

typedef struct spMatrix    /* matrix in compressed-column or triplet form */
{
    csi nzmax ;     	/* maximum number of entries */
    csi nrows ;         /* number of rows */
    csi ncols ;         /* number of columns */
    csi *cptrs ;        /* column pointers (size n+1) or col indices (size nzmax) */
    csi *rinds ;        /* row indices, size nzmax */
    csw *values ;    /* numerical values, size nzmax */
	csi nz ;        	/* # of entries in triplet matrix, -1 for compressed-col */
} spMat;

//void *spMat_calloc (csi n, size_t size) ;

//void *spMat_free (void *p) ;

void *spMat_realloc (void *p, csi n, size_t size, csi *ok) ;

spMat *spMat_spalloc (csi m, csi n, csi nzmax, csi values, csi triplet,long *dmc) ;

spMat *spMat_spfree (spMat *A,long *dmc) ;

csi spMat_sprealloc (spMat *A, csi nzmax,long *dmc) ;

//void *spMat_malloc (csi n, size_t size) ;

csi spMat_entry (spMat *T, csi r, csi c, csw v,long *dmc) ;

spMat *spMat_compress (const spMat *T,long *dmc) ;

spMat *spMat_done (spMat *C, void *w, void *x, csi ok,long *dmc) ;

double spMat_cumsum (csi *p, csi *c, csi n) ;

spMat *spMat_multiply (const spMat *A, const spMat *B,long *dmc) ;

csi spMat_scatter (const spMat *A, csi j, csw beta, csi *w, csw *x, csi mark,spMat *C, csi nz);

spMat* spMat_merge(spMat **subMatrices,csi nprocesses,long *dmc);

void spMat_difference (const spMat *A, const spMat *B, csw alpha, csw beta,csw * change);

//csi spMat_nonzero (csi i, csi j, csw aij, void *other);

csi spMat_dropzeros (spMat *A,long *dmc);

//csi spMat_fkeep (spMat *A, csi (*fkeep) (csi, csi, csw, void *), void *other);

csi spMat_fkeep (spMat *A,long *dmc);

#endif /* SPMAT_HPP_ */
