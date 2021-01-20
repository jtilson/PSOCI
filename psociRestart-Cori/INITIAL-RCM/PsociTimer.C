/**************************************************************************************

* Copyright (c) 2010 RENCI.
* All rights reserved. This program and the accompanying materials
* MAY BE available under the terms of the RENCI Open Source License
* UNC at Chapel Hill which accompanies this distribution, and is available at
* http://www.renci.org/resources/open-source-software-license

* Timing Layer Can use MPI, GA,  or Unix timer

 Classes: 

 Description:

 History:


**************************************************************************************/
/**
 *   @file PsociTimer.C
 *
 */

#include "PsociTimer.hpp"

using namespace std;

static struct timeval walltv;

double psociTime()
{
#ifdef MPITIME
  return ( MPItime() );
#elif GATIME
  return ( GAWtime() );
#else
  return ( gtod_double() );
#endif
}

double psociDiff( const double timein )
{
#ifdef MPITIME
  return( MPItime() - timein );
#elif GATIME
  return ( GAWtime() - timein );
#else
  return( gtod_double() - timein);
#endif
}


#undef MICROSEC // Thanks Gabriel for the heads up
#define MICROSEC 0.000001

double MPItime()
{
 return (MPI_Wtime());
}

#ifdef GATIME
double GAWtime()
{
 return (GA_Wtime() );
}
#endif

double gtod_double() {

  double val;
   int   ret_valt = gettimeofday( &walltv , 0);
   if (ret_valt != 0){
     int errsv = errno;
     printf("\n  error returned in gettimeofday %d \n", errsv);
     fflush(0);
   }
   val = (double) walltv.tv_sec  +  (MICROSEC * ((double) walltv.tv_usec)) ;
   return val;
}

string getTimerVersion()
{
  return (PSOCITIMER_VERSION);
}

