/* f2c.h  --  Standard Fortran to C header file */

/**  barf  [ba:rf]  2.  "He suggested using FORTRAN, and everybody barfed."

	- From The Shogakukan DICTIONARY OF NEW ENGLISH (Second edition) */

#include "slu_Cnames.h"

#ifndef F2C_INCLUDE
#define F2C_INCLUDE

#if 0
typedef long int integer; /* 64 on 64-bit machine */
typedef long int logical;
#endif

#if defined(__LAPACK_PRECISION_QUAD)
#	include <quadmath.h>
#	define M(A) A##q
	typedef __float128 quadreal;
	typedef struct { quadreal r, i; } quadcomplex;
#	define scalar __float128
#	define scalarcomplex quadcomplex
#	define dscalar __float128
#endif

typedef int integer;
typedef int logical;

typedef char *address;
typedef short int shortint;
typedef float real;
typedef double doublereal;
typedef struct { real r, i; } complex;
typedef struct { doublereal r, i; } doublecomplex;
typedef short int shortlogical;
typedef char logical1;
typedef char integer1;
/* typedef long long longint; */ /* system-dependent */

#define TRUE_ (1)
#define FALSE_ (0)

/* Extern is for use with -E */
#ifndef Extern
#define Extern extern
#endif

#define abs(x) ((x) >= 0 ? (x) : -(x))
#define dabs(x) (abs(x))
#define f2cmin(a,b) ((a) <= (b) ? (a) : (b))
#define f2cmax(a,b) ((a) >= (b) ? (a) : (b))
#define dmin(a,b) (f2cmin(a,b))
#define dmax(a,b) (f2cmax(a,b))

#if defined(__LAPACK_PRECISION_QUAD)
#	define f__cabs(r,i) qf__cabs((r),(i))
	extern scalar qf__cabs(scalar r, scalar i);
#	define pow_di(B,E) qpow_ui((B),*(E))
	extern dscalar qpow_ui(scalar *_x, integer n);
#	define myf2cmaxloc_(w,s,e,n) qf2cmaxloc_((w),*(s),*(e),n)
	extern integer qf2cmaxloc_(scalar *w, integer s, integer e, integer *n);        
#endif


#define VOID void

#endif
