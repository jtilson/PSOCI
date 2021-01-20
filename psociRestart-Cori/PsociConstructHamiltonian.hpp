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

#ifndef PSOCICONSTRICTHAMILTONIAN_H
#define PSOCICONSTRICTHAMILTONIAN_H

#include <cstdlib>
#include <cstring>

#include "PsociTimer.hpp"
#include "PsociDRAservices.hpp"
#include "PsociDRArestart.hpp"
#include "PsociGAhamiltonian.hpp"

#include "ga++.h"

using namespace std;

//Selection of H restart action
enum HSOURCE { RESTART, CONSTRUCT, NONE };

void AllocateDRA( HSOURCE state,int ndim, int maxsef, int maxsparse,  PsociDRArestart & d_cimat,PsociDRArestart & d_icicol,
		  PsociDRArestart & d_number,PsociDRArestart & d_diag_sef);

void AllocateDRA( HSOURCE state,int ndim, int maxsef, int maxsparse, int maxsparseBig, int maxwidth, PsociDRArestart & d_cimat,PsociDRArestart & d_icicol,
		  PsociDRArestart & d_number,PsociDRArestart & d_diag_sef, PsociDRArestart & d_cimat_supp, PsociDRArestart & d_icicol_supp, PsociDRArestart & d_number_supp );

void FetchHamiltonianFromDRA( HSOURCE state,  PsociDRArestart & d_cimat,PsociDRArestart & d_icicol,PsociDRArestart & d_number,PsociDRArestart & d_diag_sef, 
			      PsociGAhamiltonian & hamilton, pair<int,double> & info);

void FetchHamiltonianFromDRA( HSOURCE state,  PsociDRArestart & d_cimat,PsociDRArestart & d_icicol,PsociDRArestart & d_number,PsociDRArestart & d_diag_sef,
                               PsociDRArestart & d_cimat_supp, PsociDRArestart & d_icicol_supp, PsociDRArestart & d_number_supp, PsociGAhamiltonian & hamilton, pair<int,double> & info);

void SendHamiltonianToDRA( HSOURCE state, PsociDRArestart & d_cimat, PsociDRArestart & d_icicol, PsociDRArestart & d_number, PsociDRArestart & d_diag_sef, 
			   PsociGAhamiltonian & hamilton, pair<int,double> & info );

void SendHamiltonianToDRA( HSOURCE state, PsociDRArestart & d_cimat, PsociDRArestart & d_icicol, PsociDRArestart & d_number, PsociDRArestart & d_diag_sef, 
			    PsociDRArestart & d_cimat_supp, PsociDRArestart & d_icicol_supp, PsociDRArestart & d_number_supp, PsociGAhamiltonian & hamilton, pair<int,double> & info );

int SpecifyDRA( HSOURCE state,  PsociDRArestart & d_cimat,PsociDRArestart & d_icicol,PsociDRArestart & d_number,PsociDRArestart & d_diag_sef, int & readSparse);
int SpecifyDRA( HSOURCE state,  PsociDRArestart & d_cimat,PsociDRArestart & d_icicol,PsociDRArestart & d_number,PsociDRArestart & d_diag_sef, 
                                  PsociDRArestart & d_cimat_supp, PsociDRArestart & d_icicol_supp, PsociDRArestart & d_number_supp, int & readSparse, 
                                  int & readSparse_supp , int & readWidth );

void FetchHamiltonianFromDRA( HSOURCE state, PsociDRArestart & d_cimat, PsociDRArestart & d_icicol, PsociDRArestart & d_number, PsociDRArestart & d_diag_sef, 
			      PsociGAhamiltonian & hamilton);

void FetchHamiltonianFromDRA( HSOURCE state, PsociDRArestart & d_cimat, PsociDRArestart & d_icicol, PsociDRArestart & d_number, PsociDRArestart & d_diag_sef, 
			      PsociDRArestart & d_cimat_supp, PsociDRArestart & d_icicol_supp, PsociDRArestart & d_number_supp, PsociGAhamiltonian & hamilton);

void SendHamiltonianToDRA( HSOURCE state,  PsociDRArestart & d_cimat,PsociDRArestart & d_icicol,PsociDRArestart & d_number,PsociDRArestart & d_diag_sef, 
			   PsociGAhamiltonian & hamilton);

void SendHamiltonianToDRA( HSOURCE state,  PsociDRArestart & d_cimat,PsociDRArestart & d_icicol,PsociDRArestart & d_number,PsociDRArestart & d_diag_sef,
                           PsociDRArestart & d_cimat_supp, PsociDRArestart & d_icicol_supp, PsociDRArestart & d_number_supp, PsociGAhamiltonian & hamilton);


void InitializeDRAservice( int g_rank,  PsociDRAservices & DRAservice );

void InitializeDRAservice( int g_rank,  PsociDRAservices & DRAservice, int IOservers );

void ExitDRASERVICE(  PsociDRAservices & DRAservice );

// WARNING: The foillowiung may stil allocate MEM for H and N
void FetchICICOLFromDRA( HSOURCE state, PsociDRArestart & d_icicol,
                              PsociDRArestart & d_icicol_supp, PsociGAhamiltonian & hamilton);

void FetchICICOLFromDRA( HSOURCE state,  PsociDRArestart & d_icicol,
                               PsociDRArestart & d_icicol_supp, PsociGAhamiltonian & hamilton, pair<int,double> & info);



#endif //PSOCICONSTRICTHAMILTONIAN_H
