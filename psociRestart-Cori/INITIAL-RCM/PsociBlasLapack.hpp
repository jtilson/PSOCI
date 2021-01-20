
/**************************************************************************************

* Copyright (c) 2010,2011 RENCI.
* All rights reserved. This program and the accompanying materials
* MAY BE available under the terms of the RENCI Open Source License
* UNC at Chapel Hill which accompanies this distribution, and is available at
* http://www.renci.org/resources/open-source-software-license


 Classes: 

 Description:  Wrappers for a few of the used BLAS/LAPACK methods

 History:


**************************************************************************************/
/**
 *   @file PsociLapack.hpp
 *
 */

#ifndef PSOCIBLASLAPACK_H
#define PSOCIBLASLAPACK_H

#ifdef MKL
#include <mkl.h>
#endif

#ifdef ACML
  extern "C" {
#include <acml.h>
  }
#endif

#ifdef GOTOBLAS
  extern "C" {
      double ddot_( int * n,  double * veci, int * di, double *vecj, int * dj );
      void daxpy_( int * n, double * value, double * veci, int * di, double *vecj, int * dj );
      void dsygv_( int * itype, char * jobz, char * uplo, int * local_nsubi,  double * local_vHv, int * local_nsubj, double * local_s, int * local_nsub, double * evals, double * work, int * lwork, int & info );
}
#endif

#include <vector>

double ddotPsoci( int n, std::vector<double> & veci, int di, std::vector<double> & vecj, int dj ); 
void daxpyPsoci( int n, double value, std::vector<double> & veci, int di, std::vector<double> & vecj, int dj );
void dsygvPsoci( int itype, char jobz, char uplo, int local_nsubi,  std::vector<double> & local_vHv, int local_nsubj, std::vector<double> & local_s, int local_nsub, std::vector<double> & evals, std::vector<double> & work, int lwork, int & info );

#endif
