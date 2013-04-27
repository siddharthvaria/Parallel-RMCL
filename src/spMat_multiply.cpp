/*
 * spMat_multiply.cpp
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

/* C = A*B */

spMat* spMat_multiply (const spMat *A, const spMat *B,long *dmc) {
    csi p, j, nz = 0, anz, *Cp, *Ci, *Bp, m, n, bnz, *w, values, *Bi ;
    csw *x, *Bx, *Cx ;
    spMat *C ;
    if (!SPMAT_CSC (A) || !SPMAT_CSC (B)) return (NULL) ;      /* check inputs */
    if (A->ncols != B->nrows) return (NULL) ;
    m = A->nrows ; anz = A->cptrs [A->ncols] ;
    n = B->ncols ; Bp = B->cptrs ; Bi = B->rinds ; Bx = B->values ; bnz = Bp [n] ;
    //w = (csi*) spMat_calloc (m, sizeof (csi)) ;                    /* get workspace */
    w = (csi*) calloc(SPMAT_MAX(m,1), sizeof (csi)) ;                    /* get workspace */
    values = (A->values != NULL) && (Bx != NULL) ;
    //x = values ? (csw *) spMat_malloc (m, sizeof (csw)) : NULL ; /* get workspace */
    x = values ? (csw *) malloc (SPMAT_MAX(m,1) * sizeof (csw)) : NULL ; /* get workspace */
    C = spMat_spalloc (m, n,anz + bnz, values, 0,dmc) ;        /* allocate result */
    //cout << "Inside spMat_multiply,Dynamic memory allocated:" << *dmc << "\n";
    if (!C || !w || (values && !x)) return (spMat_done (C, w, x, 0,dmc));
    Cp = C->cptrs ;
    for (j = 0 ; j < n ; j++)
    {
        if (nz + m > C->nzmax && !spMat_sprealloc (C, 2 * (C->nzmax) + m,dmc))
        {
            return (spMat_done (C, w, x, 0,dmc)) ;             /* out of memory */
        }
        Ci = C->rinds ; Cx = C->values ;         /* C->i and C->x may be reallocated */
        Cp [j] = nz ;                   /* column j of C starts here */
        for (p = Bp [j] ; p < Bp [j+1] ; p++)
        {
            nz = spMat_scatter (A, Bi [p], Bx ? Bx [p] : 1, w, x, j+1, C, nz) ;
        }
        if (values) for (p = Cp [j] ; p < nz ; p++) Cx [p] = x [Ci [p]] ;
    }
    Cp [n] = nz ;                       /* finalize the last column of C */
    if(!spMat_sprealloc (C, 0,dmc)){                /* remove extra space from C */
    	cout << "From spMat_multiply: spMat_sprealloc failed\n";
    }
    return (spMat_done (C, w, x, 1,dmc)) ;     /* success; free workspace, return C */
}

csi spMat_scatter (const spMat *A, csi j, csw beta, csi *w, csw *x, csi mark,spMat *C, csi nz)
{
    csi i, p, *Ap, *Ai, *Ci ;
    csw *Ax ;
    if (!SPMAT_CSC (A) || !w || !SPMAT_CSC (C)) return (-1) ;     /* check inputs */
    Ap = A->cptrs ; Ai = A->rinds ; Ax = A->values ; Ci = C->rinds ;
    for (p = Ap [j] ; p < Ap [j+1] ; p++)
    {
        i = Ai [p] ;                            /* A(i,j) is nonzero */
        if (w [i] < mark)
        {
            w [i] = mark ;                      /* i is new entry in column j */
            Ci [nz++] = i ;                     /* add i to pattern of C(:,j) */
            if (x) x [i] = beta * Ax [p] ;      /* x(i) = beta*A(i,j) */
        }
        else if (x) x [i] += beta * Ax [p] ;    /* i exists in C(:,j) already */
    }
    return (nz) ;
}
