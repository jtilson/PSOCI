
/**************************************************************************************

* Copyright (c) 2010,2011 RENCI.
* All rights reserved. This program and the accompanying materials
* MAY BE available under the terms of the RENCI Open Source License
* UNC at Chapel Hill which accompanies this distribution, and is available at
* http://www.renci.org/resources/open-source-software-license


 Classes: 

 Description:  Wrappers for a few of the used BLAS/LAPACK methods
               `

 History:



**************************************************************************************/
/**
 *   @file PsociLapack.C
 *
 */

#include "PsociBlasLapack.hpp" 

using namespace std;

double ddotPsoci( int n, vector<double> & veci, int di, vector<double> & vecj, int dj ) 
{
    double value;
#ifdef MKL      
        value = ddot( &n, &veci[0], &di, &vecj[0], &dj );
#endif
#ifdef ACML
        value = ddot( n, &veci[0], di, &vecj[0], dj );
#endif
#ifdef GOTOBLAS      
        value = ddot_( &n, &veci[0], &di, &vecj[0], &dj );
#endif
#ifdef LIBSCI
        value = ddot( n, &veci[0], di, &vecj[0], dj );
#endif


    return( value );
}

void daxpyPsoci( int n, double value, vector<double> & veci, int di, vector<double> & vecj, int dj )
{
#ifdef MKL
       daxpy( &n, &value, &veci[0], &di, &vecj[0], &dj);
#endif
#ifdef ACML
       daxpy( n, value, &veci[0], di, &vecj[0], dj);
#endif
#ifdef GOTOBLAS 
       daxpy_( &n, &value, &veci[0], &di, &vecj[0], &dj);
#endif
#ifdef LIBSCI
       daxpy( n, value, &veci[0], di, &vecj[0], dj);
#endif


}


void dsygvPsoci( int itype, char jobz, char uplo, int local_nsubi, vector<double> & local_vHv, int local_nsubj, vector<double> & local_s, int local_nsub, vector<double> & evals, vector<double> & work, int lwork, int & info )
{
#ifdef MKL
      dsygv( &itype, &jobz, &uplo, &local_nsubi,  &local_vHv[0], &local_nsubj, &local_s[0], &local_nsub, &evals[0], &work[0], &lwork, &info );
#endif
#ifdef ACML
      long linfo; // Presently not using these
      dsygv( itype, jobz, uplo, local_nsubi,  &local_vHv[0], local_nsubj, &local_s[0], local_nsub, &evals[0], &info );
      //info = linfo;
#endif
#ifdef GOTOBLAS 
      int iinfo;
      dsygv_( &itype, &jobz, &uplo, &local_nsubi,  &local_vHv[0], &local_nsubj, &local_s[0], &local_nsub, &evals[0], &work[0], &lwork, iinfo );
      info = iinfo;
#endif
#ifdef LIBSCI
      long linfo; // Presently not using these
      dsygv( itype, jobz, uplo, local_nsubi,  &local_vHv[0], local_nsubj, &local_s[0], local_nsub, &evals[0], &info );
      //info = linfo;
#endif



}

