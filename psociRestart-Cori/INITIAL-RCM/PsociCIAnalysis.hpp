/**************************************************************************************

* Copyright (c) 2010,2011 RENCI.
* All rights reserved. This program and the accompanying materials
* MAY BE available under the terms of the RENCI Open Source License
* UNC at Chapel Hill which accompanies this distribution, and is available at
* http://www.renci.org/resources/open-source-software-license

* New implementation of PsociConfigs:

 Classes: 

 Description: 

 History:


**************************************************************************************/
/**
 *   @file PsociConstructHamiltonian.hpp
 *
 */

#ifndef PSOCICCIANALYSIS_H
#define PSOCICCIANALYSIS_H

#include <cstdlib>
#include <cstring>
#include <set>

#include "PsociTimer.hpp"
#include "PsociGADeterminants.hpp"
#include "PsociGAbasis.hpp"
#include "PsociGAhamiltonian.hpp"

#include "ga++.h"

using namespace std;

const int MAX_CONFIGS_PRINT = 500;

const double AUtoEV = 27.21138505; // http://physics.nist.gov/cuu/Constants/energy.html, 2011, 16 Nov

int analyzeEnergyContributions( pair<double,double> & thresh, const double fcore,  const vector<int> & l_nsef, const vector<double> & evals, const vector<double> & diags, GA::GlobalArray * g_v, PsociGAbasis & vector_set, PsociGAhamiltonian & hamiltonian, PsociGADeterminants & deters);

int energyDistributionSpatial( const double thresh, const vector<double> & deltae );

void outputSelectedRefstoCGDBG( set<string> & refs, string & excitationlevel );


bool compare_pairs( pair<int,double> left, pair<int,double> right );
bool compare_values( double left, double right );


#endif //PSOCICCIANALYSIS_H
