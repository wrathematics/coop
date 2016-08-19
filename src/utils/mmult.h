/*  Copyright (c) 2015-2016, Schmidt
    All rights reserved.
    
    Redistribution and use in source and binary forms, with or without
    modification, are permitted provided that the following conditions are met:
    
    1. Redistributions of source code must retain the above copyright notice,
    this list of conditions and the following disclaimer.
    
    2. Redistributions in binary form must reproduce the above copyright
    notice, this list of conditions and the following disclaimer in the
    documentation and/or other materials provided with the distribution.
    
    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
    "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
    TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
    PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR
    CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
    EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
    PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
    PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
    LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
    NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
    SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

#ifndef __COOP_LIB_MMULT_H__
#define __COOP_LIB_MMULT_H__


#include <stdbool.h>
#include "lapack.h"


// C ddot replica using dgemm
static inline double ddot(const int n, const double * const restrict x, const double * const restrict y)
{
  const int one = 1;
  double dot;
  
  dgemm_(&(char){'t'}, &(char){'n'}, &one, &one, &n,
    &(double){1.0}, x, &n, y, &n, &(double){0.0}, &dot, &one);
  
  return dot;
}



// upper triangle of t(x) %*% x
static inline void crossprod(const int m, const int n, const double alpha, const double * const restrict x, double *restrict c)
{
  dsyrk_(&(char){'l'}, &(char){'t'}, &n, &m, &alpha, x, &m, &(double){0.0}, c, &n);
}

static inline void tcrossprod(const int m, const int n, const double alpha, const double * const restrict x, double *restrict c)
{
  dsyrk_(&(char){'l'}, &(char){'n'}, &m, &n, &alpha, x, &m, &(double){0.0}, c, &m);
}



// dgemm wrapper
static inline void matmult(const bool transx, const bool transy, const double alpha, const int mx, const int nx, const double *const restrict x, const int my, const int ny, const double *const restrict y, double *restrict ret)
{
  // m = # rows of op(x)
  // n = # cols of op(y)
  // k = # cols of op(x)
  int im, in, ik;
  char ctransx, ctransy;
  static const double zero = 0.;
  
  ctransx = transx ? 'T' : 'N';
  ctransy = transy ? 'T' : 'N';
  
  if (transx)
  {
    im = nx;
    ik = mx;
  }
  else
  {
    im = mx;
    ik = nx;
  }
  
  in = transy ? my : ny;
  
  dgemm_(&ctransx, &ctransy, &im, &in, &ik, &alpha, x, &mx, y, &my, &zero, ret, &im);
}


#endif
