/*! \file
Copyright (c) 2003, The Regents of the University of California, through
Lawrence Berkeley National Laboratory (subject to receipt of any required 
approvals from U.S. Dept. of Energy) 

All rights reserved. 

The source code is distributed under BSD license, see the file License.txt
at the top-level directory.
*/

/*
 * -- SuperLU routine (version 5.0) --
 * Univ. of California Berkeley, Xerox Palo Alto Research Center,
 * and Lawrence Berkeley National Lab.
 * October 15, 2003
 *
 * Last update: July 10, 2015
 *
 */
#include "slu_qdefs.h"

int main(int argc, char *argv[])
{
/*
 * Purpose
 * =======
 *
 * The driver program DLINSOLX2.
 *
 * This example illustrates how to use DGSSVX to solve systems repeatedly
 * with the same sparsity pattern and similar values of matrix A.
 * In this case, the permutation vectors perm_r and perm_c are computed once.
 * The following data structures will be reused in the subsequent call to
 * DGSSVX: perm_r, perm_c, etree, L, U.
 * 
 */
    char           equed[1];
    yes_no_t       equil;
    trans_t        trans;
    SuperMatrix    A, A1, L, U;
    SuperMatrix    B, B1, X;
    NCformat       *Astore;
    NCformat       *Ustore;
    SCformat       *Lstore;
    GlobalLU_t 	   Glu;
    quadreal         *a, *a1;
    int            *asub, *xa, *asub1, *xa1;
    int            *perm_r; /* row permutations from partial pivoting */
    int            *perm_c; /* column permutation vector */
    int            *etree;
    void           *work;
    int            info, lwork, nrhs, ldx;
    int            i, j, m, n, nnz;
    quadreal         *rhsb, *rhsb1, *rhsx, *xact;
    quadreal         *R, *C;
    quadreal         *ferr, *berr;
    quadreal         u, rpg, rcond;
    mem_usage_t    mem_usage;
    superlu_options_t options;
    SuperLUStat_t stat;
    FILE 	   *fp = stdin;

    extern void    parse_command_line();

#if ( DEBUGlevel>=1 )
    CHECK_MALLOC("Enter main()");
#endif

    /* Defaults */
    lwork = 0;
    nrhs  = 1;
    equil = YES;	
    u     = 1.0;
    trans = NOTRANS;

    /* Set the default input options:
	options.Fact = DOFACT;
        options.Equil = YES;
    	options.ColPerm = COLAMD;
	options.DiagPivotThresh = 1.0;
    	options.Trans = NOTRANS;
    	options.IterRefine = NOREFINE;
    	options.SymmetricMode = NO;
    	options.PivotGrowth = NO;
    	options.ConditionNumber = NO;
    	options.PrintStat = YES;
     */
    set_default_options(&options);

    /* Can use command line input to modify the defaults. */
    parse_command_line(argc, argv, &lwork, &u, &equil, &trans);
    options.Equil = equil;
    options.DiagPivotThresh = u;
    options.Trans = trans;

    if ( lwork > 0 ) {
	work = SUPERLU_MALLOC(lwork);
	if ( !work ) {
	    ABORT("DLINSOLX: cannot allocate work[]");
	}
    }

    /* Read matrix A from a file in Harwell-Boeing format.*/
    qreadhb(fp, &m, &n, &nnz, &a, &asub, &xa);
    if ( !(a1 = quadrealMalloc(nnz)) ) ABORT("Malloc fails for a1[].");
    if ( !(asub1 = intMalloc(nnz)) ) ABORT("Malloc fails for asub1[].");
    if ( !(xa1 = intMalloc(n+1)) ) ABORT("Malloc fails for xa1[].");
    for (i = 0; i < nnz; ++i) {
        a1[i] = a[i];
	asub1[i] = asub[i];
    }
    for (i = 0; i < n+1; ++i) xa1[i] = xa[i];
    
    qCreate_CompCol_Matrix(&A, m, n, nnz, a, asub, xa, SLU_NC, SLU_Q, SLU_GE);
    Astore = A.Store;
    printf("Dimension %dx%d; # nonzeros %d\n", A.nrow, A.ncol, Astore->nnz);
    
    if ( !(rhsb = quadrealMalloc(m * nrhs)) ) ABORT("Malloc fails for rhsb[].");
    if ( !(rhsb1 = quadrealMalloc(m * nrhs)) ) ABORT("Malloc fails for rhsb1[].");
    if ( !(rhsx = quadrealMalloc(m * nrhs)) ) ABORT("Malloc fails for rhsx[].");
    qCreate_Dense_Matrix(&B, m, nrhs, rhsb, m, SLU_DN, SLU_Q, SLU_GE);
    qCreate_Dense_Matrix(&X, m, nrhs, rhsx, m, SLU_DN, SLU_Q, SLU_GE);
    xact = quadrealMalloc(n * nrhs);
    ldx = n;
    qGenXtrue(n, nrhs, xact, ldx);
    qFillRHS(trans, nrhs, xact, ldx, &A, &B);
    for (j = 0; j < nrhs; ++j)
        for (i = 0; i < m; ++i) rhsb1[i+j*m] = rhsb[i+j*m];
    
    if ( !(perm_c = intMalloc(n)) ) ABORT("Malloc fails for perm_c[].");
    if ( !(perm_r = intMalloc(m)) ) ABORT("Malloc fails for perm_r[].");
    if ( !(etree = intMalloc(n)) ) ABORT("Malloc fails for etree[].");
    if ( !(R = (quadreal *) SUPERLU_MALLOC(A.nrow * sizeof(quadreal))) )
        ABORT("SUPERLU_MALLOC fails for R[].");
    if ( !(C = (quadreal *) SUPERLU_MALLOC(A.ncol * sizeof(quadreal))) )
        ABORT("SUPERLU_MALLOC fails for C[].");
    if ( !(ferr = (quadreal *) SUPERLU_MALLOC(nrhs * sizeof(quadreal))) )
        ABORT("SUPERLU_MALLOC fails for ferr[].");
    if ( !(berr = (quadreal *) SUPERLU_MALLOC(nrhs * sizeof(quadreal))) )
        ABORT("SUPERLU_MALLOC fails for berr[].");

    /* Initialize the statistics variables. */
    StatInit(&stat);
    
    /* ------------------------------------------------------------
       WE SOLVE THE LINEAR SYSTEM FOR THE FIRST TIME: AX = B
       ------------------------------------------------------------*/
    qgssvx(&options, &A, perm_c, perm_r, etree, equed, R, C,
           &L, &U, work, lwork, &B, &X, &rpg, &rcond, ferr, berr,
           &Glu, &mem_usage, &stat, &info);

    printf("First system: dgssvx() returns info %d\n", info);

    if ( info == 0 || info == n+1 ) {

        /* This is how you could access the solution matrix. */
        quadreal *sol = (quadreal*) ((DNformat*) X.Store)->nzval;

	if ( options.PivotGrowth ) printf("Recip. pivot growth = %e\n", (double)rpg);
	if ( options.ConditionNumber )
	    printf("Recip. condition number = %e\n", (double)rcond);
        Lstore = (SCformat *) L.Store;
        Ustore = (NCformat *) U.Store;
	printf("No of nonzeros in factor L = %d\n", Lstore->nnz);
    	printf("No of nonzeros in factor U = %d\n", Ustore->nnz);
    	printf("No of nonzeros in L+U = %d\n", Lstore->nnz + Ustore->nnz - n);
    	printf("FILL ratio = %.1f\n", (float)(Lstore->nnz + Ustore->nnz - n)/nnz);

	printf("L\\U MB %.3f\ttotal MB needed %.3f\n",
	       mem_usage.for_lu/1e6, mem_usage.total_needed/1e6);
	if ( options.IterRefine ) {
            printf("Iterative Refinement:\n");
	    printf("%8s%8s%16s%16s\n", "rhs", "Steps", "FERR", "BERR");
	    for (i = 0; i < nrhs; ++i)
	      printf("%8d%8d%16e%16e\n", i+1, stat.RefineSteps, (double)ferr[i], (double)berr[i]);
	}
	fflush(stdout);

    } else if ( info > 0 && lwork == -1 ) {
        printf("** Estimated memory: %d bytes\n", info - n);
    }

    if ( options.PrintStat ) StatPrint(&stat);
    StatFree(&stat);
    Destroy_CompCol_Matrix(&A);
    Destroy_Dense_Matrix(&B);

    /* ------------------------------------------------------------
       NOW WE SOLVE ANOTHER LINEAR SYSTEM: A1*X = B1
       ONLY THE SPARSITY PATTERN OF A1 IS THE SAME AS THAT OF A.
       ------------------------------------------------------------*/
    options.Fact = SamePattern_SameRowPerm;
    StatInit(&stat); /* Initialize the statistics variables. */

    qCreate_CompCol_Matrix(&A1, m, n, nnz, a1, asub1, xa1,
                           SLU_NC, SLU_Q, SLU_GE);
    qCreate_Dense_Matrix(&B1, m, nrhs, rhsb1, m, SLU_DN, SLU_Q, SLU_GE);

    qgssvx(&options, &A1, perm_c, perm_r, etree, equed, R, C,
           &L, &U, work, lwork, &B1, &X, &rpg, &rcond, ferr, berr,
           &Glu, &mem_usage, &stat, &info);

    printf("\nSecond system: dgssvx() returns info %d\n", info);

    if ( info == 0 || info == n+1 ) {

        /* This is how you could access the solution matrix. */
        quadreal *sol = (quadreal*) ((DNformat*) X.Store)->nzval;

	if ( options.PivotGrowth ) printf("Recip. pivot growth = %e\n", (double)rpg);
	if ( options.ConditionNumber )
	    printf("Recip. condition number = %e\n", (double)rcond);
        Lstore = (SCformat *) L.Store;
        Ustore = (NCformat *) U.Store;
	printf("No of nonzeros in factor L = %d\n", Lstore->nnz);
    	printf("No of nonzeros in factor U = %d\n", Ustore->nnz);
    	printf("No of nonzeros in L+U = %d\n", Lstore->nnz + Ustore->nnz - n);
	printf("L\\U MB %.3f\ttotal MB needed %.3f\n",
	       mem_usage.for_lu/1e6, mem_usage.total_needed/1e6);
	if ( options.IterRefine ) {
            printf("Iterative Refinement:\n");
	    printf("%8s%8s%16s%16s\n", "rhs", "Steps", "FERR", "BERR");
	    for (i = 0; i < nrhs; ++i)
	      printf("%8d%8d%16e%16e\n", i+1, stat.RefineSteps, (double)ferr[i], (double)berr[i]);
	}
	fflush(stdout);
    } else if ( info > 0 && lwork == -1 ) {
        printf("** Estimated memory: %d bytes\n", info - n);
    }

    if ( options.PrintStat ) StatPrint(&stat);
    StatFree(&stat);

    SUPERLU_FREE (xact);
    SUPERLU_FREE (etree);
    SUPERLU_FREE (perm_r);
    SUPERLU_FREE (perm_c);
    SUPERLU_FREE (R);
    SUPERLU_FREE (C);
    SUPERLU_FREE (ferr);
    SUPERLU_FREE (berr);
    Destroy_CompCol_Matrix(&A1);
    Destroy_Dense_Matrix(&B1);
    Destroy_Dense_Matrix(&X);
    if ( lwork == 0 ) { /* Deallocate storage associated with L and U. */
        Destroy_SuperNode_Matrix(&L);
        Destroy_CompCol_Matrix(&U);
    } else if ( lwork > 0 ) {
        SUPERLU_FREE(work);
    }

#if ( DEBUGlevel>=1 )
    CHECK_MALLOC("Exit main()");
#endif
}

/*  
 * Parse command line options to get relaxed snode size, panel size, etc.
 */
void
parse_command_line(int argc, char *argv[], int *lwork,
                   quadreal *u, yes_no_t *equil, trans_t *trans )
{
    int c;
    extern char *optarg;

    while ( (c = getopt(argc, argv, "hl:u:e:t:")) != EOF ) {
	switch (c) {
	  case 'h':
	    printf("Options:\n");
	    printf("\t-l <int> - length of work[*] array\n");
	    printf("\t-u <int> - pivoting threshold\n");
	    printf("\t-e <0 or 1> - equilibrate or not\n");
	    printf("\t-t <0 or 1> - solve transposed system or not\n");
	    exit(1);
	    break;
	  case 'l': *lwork = atoi(optarg);
	            break;
	  case 'u': *u = atof(optarg); 
	            break;
	  case 'e': *equil = atoi(optarg); 
	            break;
	  case 't': *trans = atoi(optarg);
	            break;
  	}
    }
}
