#ifndef f77_dcomplx      /* cpp code for rename of double complex */
#define f77_dcomplx      /* fortran intrinsic functions */

#ifndef SINGLE_PRECISION
#ifndef DOUBLE_PRECISION

#ifdef NO_COPY_COMMENTS
/*  if neither SINGLE_PRECISION nor DOUBLE_PRECISION are defined, */
/*  assume DOUBLE_PRECISION  unless _CRAY is set.  */
#endif

#ifdef _CRAY
#define SINGLE_PRECISION
#else
#define DOUBLE_PRECISION
#endif

#endif  /* ifndef DOUBLE_PRECISION */
#endif  /* ifndef SINGLE_PRECISION */

#ifdef DOUBLE_PRECISION       /* workstation f77 -- double precision */

#define areal DBLE
#define AREAL DBLE
#define Areal DBLE
#define aimag DIMAG
#define AIMAG DIMAG
#define Aimag DIMAG
#define cmplx DCMPLX
#define CMPLX DCMPLX
#define Cmplx DCMPLX

#else                         /* CRAY f77 source -- single precision */

#define areal REAL
#define AREAL REAL
#define Areal REAL
#undef aimag
#undef AIMAG
#undef Aimag
#undef cmplx
#undef CMPLX
#undef Cmplx

#endif  /* DOUBLE or SINGLE_PRECISION */

#endif  /* ifndef f77_dcomplx */
