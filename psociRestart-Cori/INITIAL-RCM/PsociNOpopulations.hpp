/**************************************************************************************

* Copyright (c) 2010,2011,2012 RENCI.
* All rights reserved. This program and the accompanying materials
* MAY BE available under the terms of the RENCI Open Source License
* UNC at Chapel Hill which accompanies this distribution, and is available at
* http://www.renci.org/resources/open-source-software-license

* New implementation of PsociNOpopulations:

 Classes: 

 Description: 

 History:


**************************************************************************************/
/**
 *   @file PsociNOpopulations.hpp
 *
 */

#ifndef PSOCINOPOPULATIONS_H
#define PSOCINOPOPULATIONS_H

#include <cstdlib>
#include <cstring>
#include <string>
#include <iostream>
#include <vector>

#include "PsociGAhamiltonian.hpp"
#include "PsociProcessMOcoefs.hpp"
#include "PsociTimer.hpp"


#include "ga++.h"

//Fortran conversions

using namespace std;

bool compare_evalues( pair<double,vector<double> > left, pair<double,vector<double> > right );

const int MIN_EVAL_BOND_POPULATION = 0.05; // Do not do Mulliken bond resoltuion on orbitals with less than this occupation

enum CONDENSE { SHELL, SUBSHELL, ATOMIC  };
class PsociNOpopulations {



private:

 GA::GlobalArray * l_g_v;
 PsociProcessMOcoefs * l_mocoef;
 PsociProcessMOcoefs * l_nocoef;
 PsociGAhamiltonian * l_hamilton;

 int l_nbf;
 int l_numroots;
 int l_maxsef;
 int l_nelec;

public:

  explicit PsociNOpopulations(int numRoots, GA::GlobalArray * g_v, PsociGAhamiltonian * hamilton );

  void selectMOs( string title, int mo_unit, string mo_filename );
  void selectNOs( string title, int no_unit, int numroots, string no_filename ); // Need one file per root
  int readMOs();
  void openNOs();
  void closeNOs();
  void removeNOs();
  void removeMOs();
  long populationDensityMatrix();
  long populationDensityMatrix( pair<int,double> & info );

  long constructDensityMatrix( int numRoots, GA::GlobalArray * g_v, PsociProcessMOcoefs & coefs );
  void constructNaturalOrbitals( bool square, int nbf, int nroots, vector<vector<double> > & array, PsociProcessMOcoefs & coefs );
  void constructNaturalOrbitals( bool square, int nbf, int nroots, vector<vector<double> > & array, PsociProcessMOcoefs & coefs, pair<int,double> & info);
  void outputOccupationsSetvectors(string title,  int isym, vector<double> evals, vector<double> vectors);
  int diagonalizeNOinMObasis( int isym, int nstart, int numorbs, vector<double> & no_occupations, vector<double> & dens_vector_sym, vector<double> & array );
  void computeMOdensity( int isym, vector<double> & densityMO, vector<double> densityCI, PsociProcessMOcoefs & coefs );

  inline int twoDindex( int m, int n );
  void outputPacked2Dmatrix(bool square, int nbf, vector<double> & array, string title );
  void outputPacked1Darray(int nelem, vector<double> & array, string title );
  void outputPrintDensityAObasis( int nsym, vector<int> & nsize, vector<int> & nstart, vector<double> & density, string title );
  void replicateCIdensity(bool square, int nbf, int nroots, vector<vector<double> > & array );
  void replicateCIdensity(bool square, int nbf, int nroots, vector<vector<double> > & array, pair<int,double> & info);

  void transformNOstoAObasis( int isym, vector<double> & no_occupations, vector<double> & dens_vector_sym, vector<double> & dens_vectorAO_sym, PsociProcessMOcoefs & coefs );
  void sortNOinMObasis(int isym, vector<double> & evals, vector<double> & vectors );
  PsociProcessMOcoefs * fetchMOcoefs();
  PsociProcessMOcoefs * fetchNOcoefs();
 
  void writeNOheaders( string title );
  void fetchOrbitalFileHeader();

  void constructMullikenPopulationDensAO(int nbf, int nroots, vector<double> & densAO,  PsociProcessMOcoefs & coefs  );
  void constructMullikenPopulationNOs(int nbf, vector<double> & pk_evals, vector<double> & pk_coefs, PsociProcessMOcoefs & coefs  );
  void checkOrthonormality(int nbf, PsociProcessMOcoefs & coefs  );
  int clusterBasisLabels( CONDENSE condense, vector<string> & inlabels, vector<string> & outlabels );
  double atomicPopulations(int nstart, int nend, vector<string> & labels, vector<vector<double> > populations );
  double subshellPopulations(int nstart, int nend, vector<string> & labels, vector<vector<double> > population );
  double atomicOverlapPopulations(int nstart, int nend, vector<string> & labels, vector<vector<double> > population );

};

#endif //PSOCINOPOPULATIONS_H
