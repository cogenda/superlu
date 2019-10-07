/*! \file
Copyright (c) 2003, The Regents of the University of California, through
Lawrence Berkeley National Laboratory (subject to receipt of any required 
approvals from U.S. Dept. of Energy) 

All rights reserved. 

The source code is distributed under BSD license, see the file License.txt
at the top-level directory.
*/

/*! @file slu_qdefs.h
 * \brief Header file for real operations
 * 
 * <pre> 
 * -- SuperLU routine (version 4.1) --
 * Univ. of California Berkeley, Xerox Palo Alto Research Center,
 * and Lawrence Berkeley National Lab.
 * November, 2010
 * 
 * Global data structures used in LU factorization -
 * 
 *   nsuper: #supernodes = nsuper + 1, numbered [0, nsuper].
 *   (xsup,supno): supno[i] is the supernode no to which i belongs;
 *	xsup(s) points to the beginning of the s-th supernode.
 *	e.g.   supno 0 1 2 2 3 3 3 4 4 4 4 4   (n=12)
 *	        xsup 0 1 2 4 7 12
 *	Note: dfs will be performed on supernode rep. relative to the new 
 *	      row pivoting ordering
 *
 *   (xlsub,lsub): lsub[*] contains the compressed subscript of
 *	rectangular supernodes; xlsub[j] points to the starting
 *	location of the j-th column in lsub[*]. Note that xlsub 
 *	is indexed by column.
 *	Storage: original row subscripts
 *
 *      During the course of sparse LU factorization, we also use
 *	(xlsub,lsub) for the purpose of symmetric pruning. For each
 *	supernode {s,s+1,...,t=s+r} with first column s and last
 *	column t, the subscript set
 *		lsub[j], j=xlsub[s], .., xlsub[s+1]-1
 *	is the structure of column s (i.e. structure of this supernode).
 *	It is used for the storage of numerical values.
 *	Furthermore,
 *		lsub[j], j=xlsub[t], .., xlsub[t+1]-1
 *	is the structure of the last column t of this supernode.
 *	It is for the purpose of symmetric pruning. Therefore, the
 *	structural subscripts can be rearranged without making physical
 *	interchanges among the numerical values.
 *
 *	However, if the supernode has only one column, then we
 *	only keep one set of subscripts. For any subscript interchange
 *	performed, similar interchange must be done on the numerical
 *	values.
 *
 *	The last column structures (for pruning) will be removed
 *	after the numercial LU factorization phase.
 *
 *   (xlusup,lusup): lusup[*] contains the numerical values of the
 *	rectangular supernodes; xlusup[j] points to the starting
 *	location of the j-th column in storage vector lusup[*]
 *	Note: xlusup is indexed by column.
 *	Each rectangular supernode is stored by column-major
 *	scheme, consistent with Fortran 2-dim array storage.
 *
 *   (xusub,ucol,usub): ucol[*] stores the numerical values of
 *	U-columns outside the rectangular supernodes. The row
 *	subscript of nonzero ucol[k] is stored in usub[k].
 *	xusub[i] points to the starting location of column i in ucol.
 *	Storage: new row subscripts; that is subscripts of PA.
 * </pre>
 */
#ifndef __SUPERLU_qSP_DEFS /* allow multiple inclusions */
#define __SUPERLU_qSP_DEFS

/*
 * File name:		dsp_defs.h
 * Purpose:             Sparse matrix types and function prototypes
 * History:
 */

#ifdef _CRAY
#include <fortran.h>
#endif

/* Define my integer type int_t */
typedef int int_t; /* default */
typedef __float128 quadreal;

#include <math.h>
#include <quadmath.h>  
#include <limits.h>
#include <limits.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include "slu_Cnames.h"
#include "supermatrix.h"
#include "slu_util.h"


/* -------- Prototypes -------- */

#ifdef __cplusplus
extern "C" {
#endif

/*! \brief Driver routines */
extern void
qgssv(superlu_options_t *, SuperMatrix *, int *, int *, SuperMatrix *,
      SuperMatrix *, SuperMatrix *, SuperLUStat_t *, int *);
extern void
qgssvx(superlu_options_t *, SuperMatrix *, int *, int *, int *,
       char *, quadreal *, quadreal *, SuperMatrix *, SuperMatrix *,
       void *, int, SuperMatrix *, SuperMatrix *,
       quadreal *, quadreal *, quadreal *, quadreal *,
       GlobalLU_t *, mem_usage_t *, SuperLUStat_t *, int *);
    /* ILU */
extern void
qgsisv(superlu_options_t *, SuperMatrix *, int *, int *, SuperMatrix *,
      SuperMatrix *, SuperMatrix *, SuperLUStat_t *, int *);
extern void
qgsisx(superlu_options_t *, SuperMatrix *, int *, int *, int *,
       char *, quadreal *, quadreal *, SuperMatrix *, SuperMatrix *,
       void *, int, SuperMatrix *, SuperMatrix *, quadreal *, quadreal *,
       GlobalLU_t *, mem_usage_t *, SuperLUStat_t *, int *);


/*! \brief Supernodal LU factor related */
extern void
qCreate_CompCol_Matrix(SuperMatrix *, int, int, int, quadreal *,
		       int *, int *, Stype_t, Dtype_t, Mtype_t);
extern void
qCreate_CompRow_Matrix(SuperMatrix *, int, int, int, quadreal *,
		       int *, int *, Stype_t, Dtype_t, Mtype_t);
extern void
qCopy_CompCol_Matrix(SuperMatrix *, SuperMatrix *);
extern void
qCreate_Dense_Matrix(SuperMatrix *, int, int, quadreal *, int,
		     Stype_t, Dtype_t, Mtype_t);
extern void
qCreate_SuperNode_Matrix(SuperMatrix *, int, int, int, quadreal *, 
		         int *, int *, int *, int *, int *,
			 Stype_t, Dtype_t, Mtype_t);
extern void
qCopy_Dense_Matrix(int, int, quadreal *, int, quadreal *, int);

extern void    countnz (const int, int *, int *, int *, GlobalLU_t *);
extern void    ilu_countnz (const int, int *, int *, GlobalLU_t *);
extern void    fixupL (const int, const int *, GlobalLU_t *);

extern void    qallocateA (int, int, quadreal **, int **, int **);
extern void    qgstrf (superlu_options_t*, SuperMatrix*,
                       int, int, int*, void *, int, int *, int *, 
                       SuperMatrix *, SuperMatrix *, GlobalLU_t *,
		       SuperLUStat_t*, int *);
extern int     qsnode_dfs (const int, const int, const int *, const int *,
			     const int *, int *, int *, GlobalLU_t *);
extern int     qsnode_bmod (const int, const int, const int, quadreal *,
                              quadreal *, GlobalLU_t *, SuperLUStat_t*);
extern void    qpanel_dfs (const int, const int, const int, SuperMatrix *,
			   int *, int *, quadreal *, int *, int *, int *,
			   int *, int *, int *, int *, GlobalLU_t *);
extern void    qpanel_bmod (const int, const int, const int, const int,
                           quadreal *, quadreal *, int *, int *,
			   GlobalLU_t *, SuperLUStat_t*);
extern int     qcolumn_dfs (const int, const int, int *, int *, int *, int *,
			   int *, int *, int *, int *, int *, GlobalLU_t *);
extern int     qcolumn_bmod (const int, const int, quadreal *,
			   quadreal *, int *, int *, int,
                           GlobalLU_t *, SuperLUStat_t*);
extern int     qcopy_to_ucol (int, int, int *, int *, int *,
                              quadreal *, GlobalLU_t *);         
extern int     qpivotL (const int, const quadreal, int *, int *, 
                         int *, int *, int *, GlobalLU_t *, SuperLUStat_t*);
extern void    qpruneL (const int, const int *, const int, const int,
			  const int *, const int *, int *, GlobalLU_t *);
extern void    qreadmt (int *, int *, int *, quadreal **, int **, int **);
extern void    qGenXtrue (int, int, quadreal *, int);
extern void    qFillRHS (trans_t, int, quadreal *, int, SuperMatrix *,
			  SuperMatrix *);
extern void    qgstrs (trans_t, SuperMatrix *, SuperMatrix *, int *, int *,
                        SuperMatrix *, SuperLUStat_t*, int *);
/* ILU */
extern void    qgsitrf (superlu_options_t*, SuperMatrix*, int, int, int*,
		        void *, int, int *, int *, SuperMatrix *, SuperMatrix *,
                        GlobalLU_t *, SuperLUStat_t*, int *);
extern int     qldperm(int, int, int, int [], int [], quadreal [],
                        int [],	quadreal [], quadreal []);
extern int     ilu_qsnode_dfs (const int, const int, const int *, const int *,
			       const int *, int *, GlobalLU_t *);
extern void    ilu_qpanel_dfs (const int, const int, const int, SuperMatrix *,
			       int *, int *, quadreal *, quadreal *, int *, int *,
			       int *, int *, int *, int *, GlobalLU_t *);
extern int     ilu_qcolumn_dfs (const int, const int, int *, int *, int *,
				int *, int *, int *, int *, int *,
				GlobalLU_t *);
extern int     ilu_qcopy_to_ucol (int, int, int *, int *, int *,
                                  quadreal *, int, milu_t, quadreal, int,
                                  quadreal *, int *, GlobalLU_t *, quadreal *);
extern int     ilu_qpivotL (const int, const quadreal, int *, int *, int, int *,
			    int *, int *, int *, quadreal, milu_t,
                            quadreal, GlobalLU_t *, SuperLUStat_t*);
extern int     ilu_qdrop_row (superlu_options_t *, int, int, quadreal,
                              int, int *, quadreal *, GlobalLU_t *, 
                              quadreal *, quadreal *, int);


/*! \brief Driver related */

extern void    qgsequ (SuperMatrix *, quadreal *, quadreal *, quadreal *,
			quadreal *, quadreal *, int *);
extern void    qlaqgs (SuperMatrix *, quadreal *, quadreal *, quadreal,
                        quadreal, quadreal, char *);
extern void    qgscon (char *, SuperMatrix *, SuperMatrix *, 
		         quadreal, quadreal *, SuperLUStat_t*, int *);
extern quadreal   qPivotGrowth(int, SuperMatrix *, int *, 
                            SuperMatrix *, SuperMatrix *);
extern void    qgsrfs (trans_t, SuperMatrix *, SuperMatrix *,
                       SuperMatrix *, int *, int *, char *, quadreal *, 
                       quadreal *, SuperMatrix *, SuperMatrix *,
                       quadreal *, quadreal *, SuperLUStat_t*, int *);

extern int     sp_qtrsv (char *, char *, char *, SuperMatrix *,
			SuperMatrix *, quadreal *, SuperLUStat_t*, int *);
extern int     sp_qgemv (char *, quadreal, SuperMatrix *, quadreal *,
			int, quadreal, quadreal *, int);

extern int     sp_qgemm (char *, char *, int, int, int, quadreal,
			SuperMatrix *, quadreal *, int, quadreal, 
			quadreal *, int);
extern         quadreal qmach(char *);   /* from C99 standard, in float.h */

/*! \brief Memory-related */
extern int     qLUMemInit (fact_t, void *, int, int, int, int, int,
                            quadreal, SuperMatrix *, SuperMatrix *,
                            GlobalLU_t *, int **, quadreal **);
extern void    qSetRWork (int, int, quadreal *, quadreal **, quadreal **);
extern void    qLUWorkFree (int *, quadreal *, GlobalLU_t *);
extern int     qLUMemXpand (int, int, MemType, int *, GlobalLU_t *);

extern quadreal  *quadrealMalloc(int);
extern quadreal  *quadrealCalloc(int);
extern int     qmemory_usage(const int, const int, const int, const int);
extern int     qQuerySpace (SuperMatrix *, SuperMatrix *, mem_usage_t *);
extern int     ilu_qQuerySpace (SuperMatrix *, SuperMatrix *, mem_usage_t *);

/*! \brief Auxiliary routines */
extern void    qreadhb(FILE *, int *, int *, int *, quadreal **, int **, int **);
extern void    qreadrb(int *, int *, int *, quadreal **, int **, int **);
extern void    qreadtriple(int *, int *, int *, quadreal **, int **, int **);
extern void    qreadMM(FILE *, int *, int *, int *, quadreal **, int **, int **);
extern void    qCompRow_to_CompCol(int, int, int, quadreal*, int*, int*,
		                   quadreal **, int **, int **);
extern void    qfill (quadreal *, int, quadreal);
extern void    qinf_norm_error (int, SuperMatrix *, quadreal *);
extern quadreal  qqselect(int, quadreal *, int);


/*! \brief Routines for debugging */
extern void    qPrint_CompCol_Matrix(char *, SuperMatrix *);
extern void    qPrint_SuperNode_Matrix(char *, SuperMatrix *);
extern void    qPrint_Dense_Matrix(char *, SuperMatrix *);
extern void    qprint_lu_col(char *, int, int, int *, GlobalLU_t *);
extern int     print_quadreal_vec(char *, int, quadreal *);
extern void    qcheck_tempv(int, quadreal *);

/*! \brief BLAS */

extern int qgemm_(const char*, const char*, const int*, const int*, const int*,
                  const quadreal*, const quadreal*, const int*, const quadreal*,
		  const int*, const quadreal*, quadreal*, const int*);
extern int qtrsv_(char*, char*, char*, int*, quadreal*, int*,
                  quadreal*, int*);
extern int qtrsm_(char*, char*, char*, char*, int*, int*,
                  quadreal*, quadreal*, int*, quadreal*, int*);
extern int qgemv_(char *, int *, int *, quadreal *, quadreal *a, int *,
                  quadreal *, int *, quadreal *, quadreal *, int *);

#ifdef __cplusplus
  }
#endif

#endif /* __SUPERLU_qSP_DEFS */

