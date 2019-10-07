/*! \file
Copyright (c) 2003, The Regents of the University of California, through
Lawrence Berkeley National Laboratory (subject to receipt of any required 
approvals from U.S. Dept. of Energy) 

All rights reserved. 

The source code is distributed under BSD license, see the file License.txt
at the top-level directory.
*/

/*! @file qselect.c
 * \brief Quickselect: returns the k-th (zero-based) largest value in A[].
 *
 * <pre>
 * -- SuperLU routine (version 4.1) --
 * Lawrence Berkeley National Laboratory.
 * November, 2010
 * </pre>
 */

#define SUPERLU_MAX(x, y) 	( (x) > (y) ? (x) : (y) )
#define SUPERLU_MIN(x, y) 	( (x) < (y) ? (x) : (y) )

#include <quadmath.h>
typedef __float128 quadreal;

quadreal  qqselect(int n, quadreal A[], int k)
{
    int i, j, p;
    double val;

    k = SUPERLU_MAX(k, 0);
    k = SUPERLU_MIN(k, n - 1);
    while (n > 1)
    {
	i = 0; j = n-1;
	p = j; val = A[p];
	while (i < j)
	{
	    for (; A[i] >= val && i < p; i++);
	    if (A[i] < val) { A[p] = A[i]; p = i; }
	    for (; A[j] <= val && j > p; j--);
	    if (A[j] > val) { A[p] = A[j]; p = j; }
	}
	A[p] = val;
	if (p == k) return val;
	else if (p > k) n = p;
	else
	{
	    p++;
	    n -= p; A += p; k -= p;
	}
    }

    return A[0];
}


double dqselect(int n, double A[], int k)
{
    int i, j, p;
    double val;

    k = SUPERLU_MAX(k, 0);
    k = SUPERLU_MIN(k, n - 1);
    while (n > 1)
    {
	i = 0; j = n-1;
	p = j; val = A[p];
	while (i < j)
	{
	    for (; A[i] >= val && i < p; i++);
	    if (A[i] < val) { A[p] = A[i]; p = i; }
	    for (; A[j] <= val && j > p; j--);
	    if (A[j] > val) { A[p] = A[j]; p = j; }
	}
	A[p] = val;
	if (p == k) return val;
	else if (p > k) n = p;
	else
	{
	    p++;
	    n -= p; A += p; k -= p;
	}
    }

    return A[0];
}

float sqselect(int n, float A[], int k)
{
    int i, j, p;
    float val;

    k = SUPERLU_MAX(k, 0);
    k = SUPERLU_MIN(k, n - 1);
    while (n > 1)
    {
	i = 0; j = n-1;
	p = j; val = A[p];
	while (i < j)
	{
	    for (; A[i] >= val && i < p; i++);
	    if (A[i] < val) { A[p] = A[i]; p = i; }
	    for (; A[j] <= val && j > p; j--);
	    if (A[j] > val) { A[p] = A[j]; p = j; }
	}
	A[p] = val;
	if (p == k) return val;
	else if (p > k) n = p;
	else
	{
	    p++;
	    n -= p; A += p; k -= p;
	}
    }

    return A[0];
}
