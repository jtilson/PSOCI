/**************************************************************************************

* Copyright (c) 2010 RENCI.
* All rights reserved. This program and the accompanying materials
* MAY BE available under the terms of the RENCI Open Source License
* UNC at Chapel Hill which accompanies this distribution, and is available at
* http://www.renci.org/resources/open-source-software-license

* New implementation of PSOCI:

 Classes: 

 Description: 

 History:


**************************************************************************************/
/**
 *   @file PsociTimer.hpp
 *
 */


#ifndef PSOCITIMER_H
#define PSOCITIMER_H

#include <cstdlib>
#include <cstdio>
#include <iostream>
#include <cstring>

#include <sys/time.h>
#include <errno.h>

#include "mpi.h"

#ifdef GATIME
#include "ga.h"
double GAWtime();
#endif

const std::string PSOCITIMER_VERSION = "1.0";

double psociTime();
double psociDiff( const double timein );
double MPItime();
std::string getTimerVersion();

double gtod_double();

#endif //PSOCITIMER_H
