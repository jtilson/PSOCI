//TODO remove the superfluous GA::sync() calls used for debugging
/**************************************************************************************

* Copyright (c) 2010,2011,2012 RENCI.
* All rights reserved. This program and the accompanying materials
* MAY BE available under the terms of the RENCI Open Source License
* UNC at Chapel Hill which accompanies this distribution, and is available at
* http://www.renci.org/resources/open-source-software-license


 Classes: 

 Description: 
             
 History:

Added new preprocessing: -DCONFIGSCRAMBLE. This replace multimap with vector + scramble code


**************************************************************************************/
/**
 *   @file PsociDeterminants.C
 *
 */

#include <cmath> 
#include <set>

#include "PsociTimer.hpp"
#include "PsociDeterminants.hpp"

// Class definitions - nearly all are collectives

using namespace std;

PsociDeterminants::PsociDeterminants( PsociConfigs * configs )
{
  GA::SERVICES.sync();
  // GA::sync();
  local_configs = configs;
  local_numcfgs=0;
  local_maxdet_spatial=0;
  global_maxdet_spatial=0;
  local_maxsef_spatial=0;
  global_maxsef_spatial=0;
  local_maxdet_perspatial=0;
  global_maxdet_perspatial=0;
  local_maxsef_perspatial=0;
  global_maxsef_perspatial=0;
  printDeterminants = false;
  
  l_nelec=0;
  l_ksym=0;
  l_mxopenshells=0;
  l_nbf=0;
}

// Defaults to printing ALL LOCALLY available completed determinants after reading, assembling, phasing, and summing

void PsociDeterminants::printFinalDeterminants()
{
  printDeterminants = true;
}

//Are spins even or odd? If odd account for E symmetry when construcing dets

int PsociDeterminants::fetchSpinParity()
{
  return( local_configs->fetchSpinParity() );
}

/* Fetch the observed maximum number of open shells found in the input configs
   This can potentially be adjusted down if analysis of the input configs warrents it
*/

int PsociDeterminants::fetchMaxOpenShells()
{
  return( l_mxopenshells );
}

int PsociDeterminants::fetchNumElectrons()
{
  return( l_nelec );
}

int PsociDeterminants::fetchNumBasisFunctions()
{
  return( l_nbf );
}

int PsociDeterminants::fetchGlobalSymmetry()
{
  return( l_ksym );
}


/* New methid to support the computeConfig model of access. WE require
   external access to the object that contains the global parameters of determinants
   That is provided here.
*/
vector<JOUTFG> * PsociDeterminants::fetchJoutfg()
{
  return( &joutfg );
}

//Empty the configuration info. THis is not needed after the GA is built
/* NOTE: if using the DIRECT mode of accessing spatials must NOT destroy
   the determinant data until you are finished comstructing H
*/
void PsociDeterminants::tearDownDeterminants()
{
  if (GA::nodeid() == 0 ) cout << "Tearing down dets and configs " << endl;
  local_configs->tearDownConfigs();
  joutfg.clear();
  local_numcfgs = 0; 
}

void PsociDeterminants::printFilename()
{
  cout << "Current configuration filename is " << endl;
  local_configs->printFilename();
}

void PsociDeterminants::printParams()
{
  vector<int> params;
  local_configs->fetchParams( params );
  cout << "    Number of basis functions = " << params[0] << endl;
  cout << "          Number of electrons = " << params[1] << endl;
  cout << "              State symmetry is " << params[2] << endl;
  cout << "     Number of configurations = " << params[3] << endl;
  cout << "Maximum number of open shells = " << params[4] << endl;
  cout << "               Spin parity is = " << params[5] << endl;
  cout << endl;
}

int PsociDeterminants::maxSpatials()
{
  return( local_configs->numTotalSpatials() );
} 

int PsociDeterminants::localSpatials()
{
  return( local_numcfgs );
}


//Deprecated
int PsociDeterminants::numBasisFunctions()
{
  return( local_configs->numBasisFunctions() );
} 

//Mostly for debuggin
void PsociDeterminants::printSpatials()
{
  local_configs->printConfigurations();
}

int PsociDeterminants::fetchLocalMaxDet()
{
  return( local_maxdet_spatial );
}

int PsociDeterminants::fetchGlobalMaxDet()
{
  return( global_maxdet_spatial );
}

int PsociDeterminants::fetchLocalMaxSef()
{
  return( local_maxsef_spatial );
}

int PsociDeterminants::fetchGlobalMaxSef()
{
  return( global_maxsef_spatial );
}

int PsociDeterminants::fetchLocalMaxDetPerSpatial()
{
  return( local_maxdet_perspatial );
}

int PsociDeterminants::fetchGlobalMaxDetPerSpatial()
{
  return( global_maxdet_perspatial );
}

int PsociDeterminants::fetchLocalMaxSefPerSpatial()
{
  return( local_maxsef_perspatial );
}

int PsociDeterminants::fetchGlobalMaxSefPerSpatial()
{
  return( global_maxsef_perspatial );
}

int PsociDeterminants::fetchGlobalMaxSpatials()
{
  return( global_numcfgs );
}

void PsociDeterminants::overrideGlobalMaxSpatials( int mxspatials )
{
  global_numcfgs = mxspatials;
}

/** This method determines across the global partition who has the largest
    value for maximum det per conf ( formally known as maxdeti ). The final result
    is brdcst to all.
*/

void PsociDeterminants::assembleGlobalMaxDetSef()
{
  char op[] = "+";
  int l_maxdet = local_maxdet_spatial;
  GA::igop( &l_maxdet, 1, op );
  global_maxdet_spatial=l_maxdet;
  
  int l_maxsef = local_maxsef_spatial;
  GA::igop( &l_maxsef, 1, op );
  global_maxsef_spatial=l_maxsef;
  
  char op2[] = "max";
  int l_maxdetpS = local_maxdet_perspatial;
  GA::igop( &l_maxdetpS, 1, op2 );
  global_maxdet_perspatial=l_maxdetpS;
  
  int l_maxsefpS = local_maxsef_perspatial;
  GA::igop( &l_maxsefpS, 1, op2 );
  global_maxsef_perspatial=l_maxsefpS;

}

/*  These do not necessarily  have to be the same as returned by configs.
    They generally SHOULD but reflect the number of global spatials actually processed
    by PsociDeterminants. The calling APP must ensure all configs were processed
*/
void PsociDeterminants::assembleGlobalNumSpatials()
{
  char op[] = "+";
  //cout << "assembleGlobalNumSpatials: local_numcfgs " << local_numcfgs << endl;;
  
  int l_cfgs = local_numcfgs;
  GA::igop( &l_cfgs, 1, op );
  global_numcfgs = l_cfgs; // the number actually PROCESSED globally
}

/* Need this one because it is possible some node did not read any configs. But, they must participate
   in the creation of the GA spaces meaning they need to know a few things. */

//Somebody knows the n_elec, mxopn, etc. So simply use a max function

void PsociDeterminants::assembleGlobalJobParameters()
{
  char op[] = "max";

  int nelec = l_nelec;
  GA::igop( &nelec, 1, op );
  l_nelec = nelec;
  
  int nbf = l_nbf;
  GA::igop( &nbf, 1, op );
  l_nbf = nbf;

  int ksym = l_ksym;
  GA::igop( &ksym, 1, op );
  l_ksym = ksym;

  int mxopn = l_mxopenshells;
  GA::igop( &mxopn, 1, op );
  l_mxopenshells = mxopn;
}


//Timer wrapped version
int PsociDeterminants::fetchConfigs( int start, int stop, pair<int,double> & info )
{ 
  int status;
  int g_rank = GA::nodeid();
  double timein = psociTime();
  status = fetchConfigs( start, stop );
  info.second = psociTime() - timein;
  info.first = g_rank;
  return( status );
}


/* grab and process all configs start to stop inclusive: Using Fortran numbering: Viz., start at 1  
   Take the replicated configs, process only those from (start,stop) generating all determinents from them.
   Then ultimately push these to GA space. The calling app ensures no overlap in start/stop across cores

    This is based on a fully replicated configs object which doesn't scale
*/

//Local call
int PsociDeterminants::fetchConfigs( int start, int stop )
{
  local_maxdet_spatial = 0;
  global_maxdet_spatial = 0;
  
  vector<pair<int,string> > selcfgs;
  if ( stop < start ) {
    cerr << "fetchConfig: stop > start " << stop << " " << start << endl;
    GA::Terminate(); 
  }
#ifdef CONFIGSCRAMBLE
  GA::error("DO NOT COMPLE with -DCONFIGSCRAMBLE: possible problem in subsequent solution",-1);
  int numcfgs = local_configs->fetchConfigAlt( start, stop, selcfgs );
#else
  int numcfgs = local_configs->fetchConfig( start, stop, selcfgs );
#endif
  local_numcfgs = numcfgs;
  
#if 0
  cout << "I fetched this many " << local_numcfgs << endl;
#endif
  
  if ( local_numcfgs == 0 ) { // This is a valid answer
    return( 0 );
  }
  
  // Begin populating JOUTFG struct elements
  // NOTE: for the INDEX we need to process selcfgs in reverse ! ( probably not a big deal anyway )
  
  /* Also decrement the index as well. If not then in parallel the order/index pairs get scrambled.
     This is not a big deal until you want to write new code and need to compare fetched iconfs across
     Differently sized parallel runs
  */
 
  vector<pair<int,string> >::iterator it;
  it = selcfgs.end();
  
  if ( selcfgs.size() != (stop-start)+1 ) {
    cerr << "fetchConfig: something amiss: selcfgs not stop-start+1: Aborting: start " << start;
    cerr << " stop = " << stop << " selcfgs is " << selcfgs.size() << endl;
    GA::error("fetchConfig: something amiss: selcfgs not stop-start+1: Aborting: start ", stop);
  }
  
//Populate PsociDeterminant class variables

  l_ksym = local_configs->fetchStateSymmetry();
  l_nbf = local_configs->fetchNumBasisFtns();
  l_mxopenshells = local_configs->fetchMaxOpenShells();
  l_nelec = local_configs->fetchNumElectrons();
  
  int rev_index = stop;
  
//start and stop are absolute numbers must not allow nodes to process the same index
//The calling app must control the index selection

  if (GA::nodeid() == 0 ) cout << "Processing configurations " << endl;
  
  for (int i= start; i <= stop; ++i ) 
    {
      JOUTFG ioutfg;
      ioutfg.completed_state=0; //only triggers true when all dets/spin projections are generated
      --it;
      //      ioutfg.index = i;
      //      rev gives a cimat print identical to old code: this one must be used
      ioutfg.index = rev_index--;
      string iconf  = (*it).second;
      ioutfg.jsym = (*it).first;
      ioutfg.ksym = l_ksym;
      ioutfg.nbf = l_nbf;
      ioutfg.mx_openshells = l_mxopenshells;
      
     // Get singly and doubly occupied orbitals
     
      vector<int> singly;
      vector<int> doubly;
     
      int iopen=0;
      int iclosed=0;
      int ielec=0;
      
      for (int j=0; j<iconf.size(); ++j) {
	if (iconf[j] == '1' ) {
	  singly.push_back( j );           
	  ++iopen;
	  ++ielec;
	}
	if (iconf[j] == '2' ) {
	  doubly.push_back( j );    
	  ++iclosed;
	  ielec += 2;
	}
      }
      ioutfg.config  = iconf;
      l_nelec = ielec;
      if ( ielec != l_nelec ) {
	GA::error("ielec != nelec", ielec);
      }
      //ioutfg.nalpha = (ielec + 1)/2; //This is based on decisions made regarding the algorithm in general
      ioutfg.nelec = ielec;
      ioutfg.num_doubly = iclosed;
      ioutfg.num_singly = iopen;
      ioutfg.singly = singly;
      ioutfg.doubly = doubly;
#if USEPOW
      ioutfg.nsefi = pow(4.0, (iopen-1)/2 ) ;
      ioutfg.ndeti = pow(2.0, max( iopen-1,0 ) );
#else
      ioutfg.nsefi = (1 << (max(iopen-1,0))/2) * (1 << (max(iopen-1,0))/2);
      ioutfg.ndeti = 1 << max( iopen-1,0 );
#endif
      local_maxdet_spatial += ioutfg.ndeti; //sum local total dets 
      local_maxsef_spatial += ioutfg.nsefi; // sum local total sefs
      
      if ( ioutfg.ndeti > local_maxdet_perspatial ) local_maxdet_perspatial = ioutfg.ndeti; 
      if ( ioutfg.nsefi > local_maxsef_perspatial ) local_maxsef_perspatial = ioutfg.nsefi;

      int ntestdeti = generateDeterminants( ioutfg );

// compute phases and orbital sums (alpha/beta for each det)

     int teststatus = computePhaseAndSort( ioutfg );

#if 0
     printSpatialsData ( ioutfg );
#endif

     joutfg.push_back( ioutfg );
    } 
  return( local_numcfgs );
}

/*Special case code it grabs the configuration but DOES NOT
  update statistic about sizes, etxc.

  Are start - stop  contiguous ? 
*/
int PsociDeterminants::specialFetchConfigs( int start, int stop )
{
  if ( GA::nodeid() != 0 ) return(0);
  
  vector<pair<int,string> > selcfgs;
  int numcfgs = local_configs->fetchConfig( start, stop, selcfgs );
  
  if ( local_numcfgs == 0 ) { // This is a valid answer
    return( 0 );
  }
  local_numcfgs = numcfgs;
  //cout << GA::nodeid() <<  " specialFetchConfigs local_numcfgs is " << local_numcfgs << endl;
  
  vector<pair<int,string> >::iterator it;
  it = selcfgs.end();
  
  if ( selcfgs.size() != (stop-start)+1 ) {    
    cerr << "fetchConfig: something amiss: selcfgs not stop-start+1: Aborting: start " << start;
    cerr << " stop = " << stop << " selcfgs is " << selcfgs.size() << endl;    GA::error("fetchConfig: something amiss: selcfgs not stop-start+1: Aborting: start ", stop);
  }

  int rev_index = stop;
  
  for (int i= start; i <= stop; ++i )
    {
      JOUTFG ioutfg;
      ioutfg.completed_state=0; //only triggers true when all dets/spin projections are generated
      --it;
      ioutfg.index = rev_index--;
      string iconf  = (*it).second;
      ioutfg.jsym = (*it).first;
      ioutfg.ksym = l_ksym;
      ioutfg.nbf = l_nbf;
      ioutfg.mx_openshells = l_mxopenshells;
      
      vector<int> singly;
      vector<int> doubly;
      
      int iopen=0;
      int iclosed=0;
      int ielec=0;
      
      for (int j=0; j<iconf.size(); ++j) {
	if (iconf[j] == '1' ) {
	  singly.push_back( j );           
	  ++iopen;
	  ++ielec;
	}
	if (iconf[j] == '2' ) {
	  doubly.push_back( j );    
	  ++iclosed;
	  ielec += 2;
	}
      }
      ioutfg.config  = iconf;
      l_nelec = ielec;
      if ( ielec != l_nelec ) {
	GA::error("ielec != nelec", ielec);
      }
      ioutfg.nelec = ielec;
      ioutfg.num_doubly = iclosed;
      ioutfg.num_singly = iopen;
      ioutfg.singly = singly;
      ioutfg.doubly = doubly;
#if USEPOW
      ioutfg.nsefi = pow(4.0, (iopen-1)/2 ) ;
      ioutfg.ndeti = pow(2.0, max( iopen-1,0 ) );
#else
      ioutfg.nsefi = (1 << (max(iopen-1,0))/2) * (1 << (max(iopen-1,0))/2);
      ioutfg.ndeti = 1 << max( iopen-1,0 );
#endif
      int ntestdeti = generateDeterminants( ioutfg );
      int teststatus = computePhaseAndSort( ioutfg );
      
      joutfg.push_back( ioutfg );
    }
  return( local_numcfgs );
}


/*Special case code it grabs the configuration but DOES NOT
  update statistic about sizes, etc.

  THis is a complicated approach. We can't replicate configs to all nodes for big jobs.
  So we have node=0, carve up configs into contiguous blocks ( or empty ones). These get
  passed ot a remote node and used to initialize deters.

  We can then simply set start=1, stop =length, BUT to get the absolute index
  correct ro GA, we use abs_index to initialize rev_index;

  in this case we do not need start/stop since we always grab all local configs
*/

int PsociDeterminants::specialDistributedFetchConfigs()
{
//Populate PsociDeterminant class variables

  l_ksym = local_configs->fetchStateSymmetry();
  l_nbf = local_configs->fetchNumBasisFtns();
  l_mxopenshells = local_configs->fetchMaxOpenShells();
  l_nelec = local_configs->fetchNumElectrons();

  assembleGlobalJobParameters();
  
  //  if ( GA::nodeid() == 0 ) cout << "specialDistributedFetchConfigs "<< l_nbf <<" "<<l_mxopenshells<<" "<<l_nelec<<endl;
  vector<pair<int,string> > selcfgs;
  vector<int> indexes;

  /* Need a new version that also grabs the absolute indexes and LOCAL
     data objects
  */
  
  int numcfgs = local_configs->fetchConfigGA( selcfgs, indexes ); 
 
  if ( numcfgs == 0 ) return(0);
  
  local_numcfgs = numcfgs; 
  vector<pair<int,string> >::iterator it;
  
  it = selcfgs.begin();
  int stop =  selcfgs.size();

//No need to do this in reverse order. We already have the indexes
  
  int rev_index = 0;
  
  for (int i= 0; i < stop; ++i )
    {
      //rindex = --rev_index;
      JOUTFG ioutfg;
      ioutfg.completed_state=0; //only triggers true when all dets/spin projections are generated
      //--it;
      
      ioutfg.index = indexes[ i ]; // I expect them to be contiguous but ....
      
      string iconf  = (*it).second;
      
      ioutfg.jsym = (*it).first;
      ioutfg.ksym = l_ksym;
      ioutfg.nbf = l_nbf;
      ioutfg.mx_openshells = l_mxopenshells;
      
      vector<int> singly;
      vector<int> doubly;
      
      int iopen=0;
      int iclosed=0;
      int ielec=0;
      
      for (int j=0; j<iconf.size(); ++j) {
	if (iconf[j] == '1' ) {
	  singly.push_back( j );           
	  ++iopen;
	  ++ielec;
	}
	if (iconf[j] == '2' ) {
	  doubly.push_back( j );    
	  ++iclosed;
	  ielec += 2;
	}
      }
      ioutfg.config  = iconf;
      l_nelec = ielec;
      if ( ielec != l_nelec ) {
	GA::error("ielec != nelec", ielec);
      }
      ioutfg.nelec = ielec;
      ioutfg.num_doubly = iclosed;
      ioutfg.num_singly = iopen;
      ioutfg.singly = singly;
      ioutfg.doubly = doubly;
#if USEPOW
      ioutfg.nsefi = pow(4.0, (iopen-1)/2 ) ;
      ioutfg.ndeti = pow(2.0, max( iopen-1,0 ) );
#else
      ioutfg.nsefi = (1 << (max(iopen-1,0))/2) * (1 << (max(iopen-1,0))/2);
      ioutfg.ndeti = 1 << max( iopen-1,0 );
#endif
      int ntestdeti = generateDeterminants( ioutfg );
      int teststatus = computePhaseAndSort( ioutfg );
      joutfg.push_back( ioutfg );
      ++it;
    }
  return( local_numcfgs );
}

/* slight modification to allow populating local phase data structure 
   a this stage we presume replicated_phases is already sized ( l_maxSpatials ) 
   and set to zero. It wil be indexed as (joutfg.index)-1 = value
*/
int PsociDeterminants::specialDistributedFetchConfigs( vector<char> & replicated_phases )
{
//Populate PsociDeterminant class variables

  l_ksym = local_configs->fetchStateSymmetry();
  l_nbf = local_configs->fetchNumBasisFtns();
  l_mxopenshells = local_configs->fetchMaxOpenShells();
  l_nelec = local_configs->fetchNumElectrons();

  assembleGlobalJobParameters();
  
  //  if ( GA::nodeid() == 0 ) cout << "specialDistributedFetchConfigs "<< l_nbf <<" "<<l_mxopenshells<<" "<<l_nelec<<endl;
  vector<pair<int,string> > selcfgs;
  vector<int> indexes;

  /* Need a new version that also grabs the absolute indexes and LOCAL
     data objects
  */
  
  int numcfgs = local_configs->fetchConfigGA( selcfgs, indexes ); 
 
  if ( numcfgs == 0 ) return(0);
  
  local_numcfgs = numcfgs; 
  vector<pair<int,string> >::iterator it;
  
  it = selcfgs.begin();
  int stop =  selcfgs.size();

//No need to do this in reverse order. We already have the indexes
  
  int rev_index = 0;
  
  if ( replicated_phases.size() != maxSpatials() ) {
     cout<<" wrong size for replicated_phases" << replicated_phases.size()<<" "<<maxSpatials() << endl;
     GA::error(" wrong size for replicated_phases", replicated_phases.size() );
  }

  for (int i= 0; i < stop; ++i )
    {
      //rindex = --rev_index;
      JOUTFG ioutfg;
      ioutfg.completed_state=0; //only triggers true when all dets/spin projections are generated
      //--it;
      
      ioutfg.index = indexes[ i ]; // I expect them to be contiguous but ....
      
      string iconf  = (*it).second;
      
      ioutfg.jsym = (*it).first;
      ioutfg.ksym = l_ksym;
      ioutfg.nbf = l_nbf;
      ioutfg.mx_openshells = l_mxopenshells;
      
      vector<int> singly;
      vector<int> doubly;
      
      int iopen=0;
      int iclosed=0;
      int ielec=0;
      
      for (int j=0; j<iconf.size(); ++j) {
	if (iconf[j] == '1' ) {
	  singly.push_back( j );           
	  ++iopen;
	  ++ielec;
	}
	if (iconf[j] == '2' ) {
	  doubly.push_back( j );    
	  ++iclosed;
	  ielec += 2;
	}
      }
      ioutfg.config  = iconf;
      l_nelec = ielec;
      if ( ielec != l_nelec ) {
	GA::error("ielec != nelec", ielec);
      }
      ioutfg.nelec = ielec;
      ioutfg.num_doubly = iclosed;
      ioutfg.num_singly = iopen;
      ioutfg.singly = singly;
      ioutfg.doubly = doubly;
#if USEPOW
      ioutfg.nsefi = pow(4.0, (iopen-1)/2 ) ;
      ioutfg.ndeti = pow(2.0, max( iopen-1,0 ) );
#else
      ioutfg.nsefi = (1 << (max(iopen-1,0))/2) * (1 << (max(iopen-1,0))/2);
      ioutfg.ndeti = 1 << max( iopen-1,0 );
#endif
      int ntestdeti = generateDeterminants( ioutfg );
      int teststatus = computePhaseAndSort( ioutfg );
      joutfg.push_back( ioutfg );
      ++it;

#ifndef REPLICATEDETVALUE 
    GA::error("PsociDeterminants::specialDistributedFetchConfigs(replicated)  should only be used with -DREPLICATEDETVALUE",1);
#endif
// phases are always the same for all contributing detemrinantes AND either 0 or 1.
#ifdef REPLICATEDETVALUE

  //  replicated_phases[ ioutfg.index-1 ] = ioutfg.phases[0]; // there is always one or more of them 
  //  int test = ioutfg.spin_projection[0];
  //  replicated_phases[ ioutfg.index-1 ] = (test % 2);
     replicated_phases[ ioutfg.index-1 ] = (char) iopen;

    //phases do not seem connected to projections SO use even/odd indicators instead
#endif
    }
  return( local_numcfgs );
}

/* override valuue for local_numcfgs
*/
void PsociDeterminants::setLocalNumcfgs( int numcfgs )
{
  local_numcfgs = numcfgs;
}

/* Method to basically get the set of nsef.ndeti for a config list
   Not the most efficient way but it works and is not a big deal.
   
   A parallel method
*/
int PsociDeterminants::findParams()
{
  local_maxdet_spatial = 0;
  global_maxdet_spatial = 0;
  local_maxsef_spatial = 0;
  global_maxsef_spatial = 0;
  
  int stop = local_configs->localVecConfigsSize(); //how many do I have ?
  if ( stop < 1 ) return(0);
  int start = 1;
  
  vector<pair<int,string> > selcfgs;
  
  if ( stop < start ) {
    cerr << "findParams: stop > start " << stop << " " << start << endl;
    GA::Terminate(); 
  }
  
  local_numcfgs = 0;
  
  selcfgs.clear(); 
  int numcfgs = local_configs->fetchConfigAlt( start, stop, selcfgs );
  local_numcfgs += numcfgs;
  
  l_ksym = local_configs->fetchStateSymmetry();
  l_nbf = local_configs->fetchNumBasisFtns();
  l_mxopenshells = local_configs->fetchMaxOpenShells();
  l_nelec = local_configs->fetchNumElectrons();

    for (int i= 0; i < stop; ++i )
      {
	JOUTFG ioutfg;
	ioutfg.completed_state=0; //only triggers true when all dets/spin projections are generated
	ioutfg.index = i;
	//      rev gives a cimat print identical to old code: this one must be used
	ioutfg.index = i; // dummy
	
	//TODO fix me this is ridiculous
	
	string iconf  = selcfgs[i].second;
	ioutfg.jsym = selcfgs[i].first;
	ioutfg.ksym = l_ksym;
	ioutfg.nbf = l_nbf;
	ioutfg.mx_openshells = l_mxopenshells;
	
	// Get singly and doubly occupied orbitals
	
	vector<int> singly;
	vector<int> doubly;
	
	int iopen=0;
	int iclosed=0;
	int ielec=0;
	
	for (int j=0; j<iconf.size(); ++j) {
	  if (iconf[j] == '1' ) {
	    singly.push_back( j );           
	    ++iopen;
	    ++ielec;
	  }
	  if (iconf[j] == '2' ) {
	    doubly.push_back( j );    
	    ++iclosed;
	    ielec += 2;
	  }
	}
	ioutfg.config  = iconf;
	l_nelec = ielec;
	if ( ielec != l_nelec ) {
	  GA::error("ielec != nelec", ielec);
	}
	//ioutfg.nalpha = (ielec + 1)/2; //This is based on decisions made regarding the algorithm in general
	ioutfg.nelec = ielec;
	ioutfg.num_doubly = iclosed;
	ioutfg.num_singly = iopen;
	//ioutfg.singly = singly;
	//ioutfg.doubly = doubly;
#if USEPOW
	ioutfg.nsefi = pow(4.0, (iopen-1)/2 ) ;
	ioutfg.ndeti = pow(2.0, max( iopen-1,0 ) );
#else
	ioutfg.nsefi = (1 << (max(iopen-1,0))/2) * (1 << (max(iopen-1,0))/2);
	ioutfg.ndeti = 1 << max( iopen-1,0 );
#endif
	local_maxdet_spatial += ioutfg.ndeti; //sum local total dets 
	local_maxsef_spatial += ioutfg.nsefi; // sum local total sefs
	
	if ( ioutfg.ndeti > local_maxdet_perspatial ) local_maxdet_perspatial = ioutfg.ndeti; 
	if ( ioutfg.nsefi > local_maxsef_perspatial ) local_maxsef_perspatial = ioutfg.nsefi;
      } 
    return( local_numcfgs );
}

/* This is the entry point to return the set of determinants for the
   indicted range of configurations.
*/
/* NOTE: configs MUST be pre replicated on all cores. This is an alternative to the "fetch"
   method. Instead of fetching determinant data from GA, we simply recompute all determinants
   for the list of configs (start,stop) inclusive. Data are returned in the vector of JOUTFG objects.
   
   This is a tricky method to use as it requires thr alternative constructor and generate method
   in PsociGADeterminents.
   
   Lastly, we use too much memory since the GA space of determinents continues to exist for use by 
   downstream analysis methods. 
*/

int PsociDeterminants::computeConfigs( int start, int stop, vector<JOUTFG> & toutfg )
{
  vector<pair<int,string> > selcfgs;

#ifdef DETAILEDCHECK
  if ( stop < start ) {
    cerr << "fetchConfig: stop > start " << stop << " " << start << endl;
    GA::Terminate();
  }
#endif
//int numcfgs = local_configs->fetchConfig( start, stop, selcfgs ); // Slower
  int numcfgs = local_configs->fetchConfigAlt( start, stop, selcfgs );
  local_numcfgs = numcfgs;
  
#if 0
  cout << "computeConfig: I fetched this many into selcfgs " << local_numcfgs << endl;
#endif
  
// Begin populating JOUTFG struc elements
// NOTE: for the INDEX we need to process selcfgs in reverse ! ( probaly not a big deal anyway )
/* In contrast to the fetch2D methods here we DO NOT want to proceed in reverse order
   When provided a list on indexes simply grab them in order.
*/

#ifdef DETAILEDCHECK
  if ( selcfgs.size() != (stop-start)+1 ) {
    cerr << "computerConfig: something amiss: selcfgs not stop-start+1: Aborting: start " << start;
    cerr << " stop = " << stop << " selcfgs is " << selcfgs.size() << endl;
    GA::error("computeConfig: something amiss: selcfgs not stop-start+1: Aborting: start ", stop);
  }
#endif

//Populate PsociDeterminant class variables

  l_ksym = local_configs->fetchStateSymmetry();
  l_nbf = local_configs->fetchNumBasisFtns();
  l_mxopenshells = local_configs->fetchMaxOpenShells();
  l_nelec = local_configs->fetchNumElectrons();
  
//start and stop are absolute numbers must not allow nodes to process the same index
//The calling app must control the index selection
//if (GA::nodeid() == 0 ) cout << "computeConfig: Processing configurations " << endl;
/* NOTE we need to transpose the order coming back from fetchConfigs
*/


//  int index = stop - start; // swap the order in a relative way
//  int index=selcfgs.size();

  int index = -1;
  
  for (int i= start; i <= stop; ++i )
    {
      JOUTFG ioutfg; //Note phase method APPENDS data to this
      
      ++index;
      ioutfg.completed_state=0; //only triggers true when all dets/spin projections are generted
      ioutfg.index = i; // Need the real value here
      
      string iconf  = selcfgs[index].second;
      ioutfg.jsym = selcfgs[index].first;
      ioutfg.ksym = local_configs->fetchStateSymmetry();
      ioutfg.nbf = local_configs->fetchNumBasisFtns();
      ioutfg.mx_openshells = local_configs->fetchMaxOpenShells();
      
      // Get singly and doubly occupied orbitals
      
      vector<int> singly;
      vector<int> doubly;
      
      int iopen=0;
      int iclosed=0;
      int ielec=0;
      
      for (int j=0; j<iconf.size(); ++j) {
        if (iconf[j] == '1' ) {
          singly.push_back( j );
          ++iopen;
          ++ielec;
        }
        if (iconf[j] == '2' ) {
          doubly.push_back( j );
          ++iclosed;
          ielec += 2;
        }
      }
      
      ioutfg.config  = iconf;
      //ioutfg.nalpha = (ielec + 1)/2; //This is based on decisions made regarding the algorithm in general
      ioutfg.nelec = ielec;
      ioutfg.num_doubly = iclosed;
      ioutfg.num_singly = iopen;
      ioutfg.singly = singly;
      ioutfg.doubly = doubly;
#if USEPOW
      ioutfg.nsefi = pow(4.0, (iopen-1)/2 ) ;
      ioutfg.ndeti = pow(2.0, max( iopen-1,0 ) );
#else
      ioutfg.nsefi = (1 << (max(iopen-1,0))/2) * (1 << (max(iopen-1,0))/2);
      ioutfg.ndeti = 1 << max( iopen-1,0 );
#endif
      
      // computePhaseAndSort "pushes" data into the followig arrays
      
      int ntestdeti = generateDeterminants( ioutfg );
      int teststatus = computePhaseAndSort( ioutfg );
#if 0
  printSpatialsData ( ioutfg );
#endif
  
  toutfg.push_back( ioutfg );
    }
  return( local_numcfgs );
}

/* This is the entry point to return the ONE determinant
*/
/* NOTE: configs MUST be pre replicated on all cores. THis is the simpler version of getting determinants. In this
   version we assume only a SINGLE configuration worth of data are required. It is returned in toutfg.
*/

int PsociDeterminants::computeConfigs( int index, JOUTFG & toutfg )
{
  double temp1;
  temp1 = psociTime();
  vector<pair<int,string> > selcfgs;
//  int numcfgs = local_configs->fetchConfig( index, index, selcfgs ); // SLower
  int numcfgs = local_configs->fetchConfigAlt( index, index, selcfgs );
  local_numcfgs = numcfgs;
  cout << " get configs " << psociTime() - temp1<<endl;
  
#if 0
  cout << "computeConfig: I fetched this many " << local_numcfgs << endl;
#endif

#ifdef DETAILEDCHECK
  if ( selcfgs.size() != 1 ) {
    cerr << "compuyeConfig: Single det: something amiss: selcfgs not stop-start+1: Aborting: index " << index;
    GA::error("computeConfig: Single det: something amiss: selcfgs not index: Aborting: ", index);
  }
#endif

//Populate PsociDeterminant class variables
/*
  l_ksym = local_configs->fetchStateSymmetry();
  l_nbf = local_configs->fetchNumBasisFtns();
  l_mxopenshells = local_configs->fetchMaxOpenShells();
  l_nelec = local_configs->fetchNumElectrons();
  int rev_index = index;
*/

  //start and stop are absolute numbers must not allow nodes to process the same index
  //The calling app must control the index selection

//  if (GA::nodeid() == 0 ) cout << "computeConfig: Processing configurations " << endl;

  toutfg.completed_state=0; //only triggers true when all dets/spin projections are generted
  toutfg.index = index;
  string iconf  = selcfgs[0].second;
  toutfg.jsym = selcfgs[0].first;
  toutfg.ksym = local_configs->fetchStateSymmetry();
  toutfg.nbf = local_configs->fetchNumBasisFtns();
  toutfg.mx_openshells = local_configs->fetchMaxOpenShells();

  // Get singly and doubly occupied orbitals
  
  vector<int> singly;
  vector<int> doubly;
  
  int iopen=0;
  int iclosed=0;
  int ielec=0;
  
  for (int j=0; j<iconf.size(); ++j) {
    if (iconf[j] == '1' ) {
      singly.push_back( j );
      ++iopen;
      ++ielec;
    }
    if (iconf[j] == '2' ) {
      doubly.push_back( j );
      ++iclosed;
      ielec += 2;
    }
  }
  toutfg.config  = iconf;
  toutfg.nelec = ielec;
  toutfg.num_doubly = iclosed;
  toutfg.num_singly = iopen;
  toutfg.singly = singly;
  toutfg.doubly = doubly;
#if USEPOW
  toutfg.nsefi = pow(4.0, (iopen-1)/2 ) ;
  toutfg.ndeti = pow(2.0, max( iopen-1,0 ) );
#else
  toutfg.nsefi = (1 << (max(iopen-1,0))/2) * (1 << (max(iopen-1,0))/2);
  toutfg.ndeti = 1 << max( iopen-1,0 );
#endif

// THis is tiny
  temp1 = psociTime();
  int ntestdeti = generateDeterminants( toutfg );
  //cout << "GENERATE time " << psociTime()-temp1<<endl;
  
  // compute phases and orbital sums (alpha/beta for each det)
  
// computePhaeAndSort "pushes" data into the followig arrays

  toutfg.phases.clear();
  toutfg.sumsa.clear();
  toutfg.sumsb.clear();

  temp1 = psociTime();
  int teststatus = computePhaseAndSort( toutfg );
  //cout << "PHASE time " << psociTime()-temp1<<endl;

#if 0
  printSpatialsData ( toutfg );
#endif
  return( local_numcfgs );
}

/* A new method based on the orbmap concept. Here we simply pass in a vector
   containing the orbmap and construct the determinates for it

   for this routine sefcfgs is not required nor is internal conversion to vector
*/
int PsociDeterminants::computeConfigs( vector<COMPRESS_SIZE> & orbmap, JOUTFG & toutfg )
{

    double temp1;
//  temp1 = psociTime();
//  vector<pair<int,string> > selcfgs;
//  int numcfgs = local_configs->fetchConfig( index, index, selcfgs ); // SLower
//  int numcfgs = local_configs->fetchConfigAlt( index, index, selcfgs );
//  local_numcfgs = numcfgs;
//  cout << " get configs " << psociTime() - temp1<<endl;
  
#if 0
  cout << "computeConfig: I fetched this many " << local_numcfgs << endl;
#endif

#ifdef DETAILEDCHECK
//  if ( selcfgs.size() != 1 ) {
//    cerr << "compuyeConfig: Single det: something amiss: selcfgs not stop-start+1: Aborting: index " << index;
//    GA::error("computeConfig: Single det: something amiss: selcfgs not index: Aborting: ", index);
//  }
#endif

//Populate PsociDeterminant class variables
/*
  l_ksym = local_configs->fetchStateSymmetry();
  l_nbf = local_configs->fetchNumBasisFtns();
  l_mxopenshells = local_configs->fetchMaxOpenShells();
  l_nelec = local_configs->fetchNumElectrons();
  int rev_index = index;
*/

  //start and stop are absolute numbers must not allow nodes to process the same index
  //The calling app must control the index selection

#if 0
  if (GA::nodeid() == 0 ) cout << "computeConfig: Processing configurations : index is" << orbmap[l_nbf+1] << endl;
#endif

  toutfg.completed_state=0; //only triggers true when all dets/spin projections are generted
  toutfg.jsym = orbmap[l_nbf];
  toutfg.index = orbmap[l_nbf+1];
  toutfg.ksym = local_configs->fetchStateSymmetry();
  toutfg.nbf = local_configs->fetchNumBasisFtns();
  toutfg.mx_openshells = local_configs->fetchMaxOpenShells();

  // Get singly and doubly occupied orbitals
  
  vector<int> singly;
  vector<int> doubly;
  
  int iopen=0;
  int iclosed=0;
  int ielec=0;
  
  for (int j=0; j<l_nbf; ++j) {
    if (orbmap[j] == 1 ) {
      singly.push_back( j );
      ++iopen;
      ++ielec;
    }
    if (orbmap[j] == 2 ) {
      doubly.push_back( j );
      ++iclosed;
      ielec += 2;
    }
  }
//  toutfg.config  = iconf;
  toutfg.nelec = ielec;
  toutfg.num_doubly = iclosed;
  toutfg.num_singly = iopen;
  toutfg.singly = singly;
  toutfg.doubly = doubly;
#if USEPOW
  toutfg.nsefi = pow(4.0, (iopen-1)/2 ) ;
  toutfg.ndeti = pow(2.0, max( iopen-1,0 ) );
#else
  toutfg.nsefi = (1 << (max(iopen-1,0))/2) * (1 << (max(iopen-1,0))/2);
  toutfg.ndeti = 1 << max( iopen-1,0 );
#endif

// THis is tiny
  temp1 = psociTime();
  int ntestdeti = generateDeterminants( toutfg );
  //cout << "GENERATE time " << psociTime()-temp1<<endl;
  
  // compute phases and orbital sums (alpha/beta for each det)
  
// computePhaeAndSort "pushes" data into the followig arrays

  toutfg.phases.clear();
  toutfg.sumsa.clear();
  toutfg.sumsb.clear();

  temp1 = psociTime();
  int teststatus = computePhaseAndSort( toutfg );
  //cout << "PHASE time " << psociTime()-temp1<<endl;

// resize off the last entries of orbmap and construct a stringified version
  orbmap.resize( l_nbf );

  string str_config;
  local_configs->convertIntArrayToString( orbmap, str_config );
  toutfg.config = str_config;
  
#if 0
  printSpatialsData ( toutfg );
#endif
  return( local_numcfgs );
}

/* A useful debugging method for grabbing one of the existing configs and
   dumping out the specifics

   This method is a pain since the index that requests a print must start at 1 and
   end end at the number of local configs. 
*/

void PsociDeterminants::printSpatialsData( int iconf )
{
  if ( iconf > local_numcfgs ) {
    cerr << " printSpatialsData: Warning: out of bounds iconfig number. Only have " << local_numcfgs ;
    cerr << " Information for index number " << iconf << " was provided: Skipping " << endl;
    return;
  }
  
  if ( joutfg.size() != local_numcfgs ) {
    cerr << "printSpatialsData: State inconsistency local nums reports" << local_numcfgs ;
    cerr << " but joutfg list ids of size " << joutfg.size() << endl;
    GA::Terminate();
  }
  
  JOUTFG ioutfg = joutfg[ iconf  - 1 ];
  printSpatialsData ( ioutfg );
}

/* I don't printout the local variables ( e.g.,l_nbf) because we could
   call this after a teardown and the local values would not be set
   leading to confusing printouts
*/

void PsociDeterminants::printSpatialsData( JOUTFG & ioutfg )
{
  int g_rank = GA::nodeid();
  if ( g_rank == 0 ) {
  cout << " completed  = " << ioutfg.completed_state << endl;
  cout << " index      = " << ioutfg.index << endl;
  cout << " nbf        = " << ioutfg.nbf << endl;
  cout << " jsym       = " << ioutfg.jsym << endl;
  cout << " ksym       = " << ioutfg.ksym << endl;
  cout << " mxopn      = " << ioutfg.mx_openshells << endl;
  cout << " nelec      = " << ioutfg.nelec << endl;
  //cout << " nalpha     = " << ioutfg.nalpha << endl;
  cout << " num_doubly = " << ioutfg.num_doubly << endl;
  cout << " num_singly = " << ioutfg.num_singly << endl;
  cout << " config: " << ioutfg.config << endl;
  cout << " nsefi      = " << ioutfg.nsefi << endl;
  cout << " ndeti      = " << ioutfg.ndeti << endl;
  }

  vector<int>::iterator it;
  
  cout << "All reported orbital indexes have been shifted: A one was added. ";
  cout << " Internally, however, all orbitals begin at zero " << endl;
  
  cout << "List of doubly occupied orbs: " << endl;
  for ( it = ioutfg.doubly.begin(); it != ioutfg.doubly.end(); ++it ) {
    cout << " " << (*it)+1;
  }
  cout << endl << "List of singly occupied orbs: " << endl;
  for ( it = ioutfg.singly.begin(); it != ioutfg.singly.end(); ++it ) {
    cout << " " << (*it)+1;
  }
  cout << endl;
  
  if ( ioutfg.alpha_spins.size() != 0 ) {
    int ndet=0;
    vector<int>::iterator sit;
    sit = ioutfg.spin_projection.begin();
    cout << "List of alpha orbitals for this determinant " << endl;
    vector<vector<int> >::iterator vit;
    vector<int>::iterator inner;
    for (vit = ioutfg.alpha_spins.begin(); vit != ioutfg.alpha_spins.end(); ++vit) {
      ++ndet;
      cout << "Determinant number: " << ndet <<" spin projection is " << (*sit) << " Number of alpha orbitals: " << (*vit).size()  <<  endl;
      for (inner=(*vit).begin(); inner !=(*vit).end(); ++inner){
	cout << " " << (*inner)+1;
      }   
      cout << endl;
      ++sit;
    } 
  } else {
    cout << "No alpha_spin orbitals have been specified " << endl;
  }

  if ( ioutfg.beta_spins.size() != 0 ) {
    int ndet=0;
    vector<int>::iterator sit;
    sit = ioutfg.spin_projection.begin();
    cout << "List of beta orbitals for this determinant " << endl;
    vector<vector<int> >::iterator vit;
    vector<int>::iterator inner;
    for (vit = ioutfg.beta_spins.begin(); vit != ioutfg.beta_spins.end(); ++vit) {
      ++ndet;
      cout << "Determinant number: " << ndet <<" spin projection is " << (*sit)  << " Number of beta orbitals: " << (*vit).size()  <<  endl;
      for (inner=(*vit).begin(); inner !=(*vit).end(); ++inner){
	cout << " " << (*inner)+1;
      } 
      cout << endl;
      ++sit;
    } 
  } else {
    cout << "No beta_spin orbitals have been specified " << endl;
  }
  
  if ( ioutfg.phases.size() != 0 ) {
    cout << "Determinant phases are " << endl;
    for (it=(ioutfg.phases).begin();  it !=ioutfg.phases.end(); ++it){
      cout << " " << (*it);
    }
    cout << endl;
  } else {
    cout << "No phases have been specified " << endl;
  }
  
   if ( ioutfg.sumsa.size() != 0 ) {
     //vector<int>::iterator it;
     cout << "Alpha orbital sums are " << endl;
     for (it=(ioutfg.sumsa).begin();  it !=ioutfg.sumsa.end(); ++it){
       cout << " " << (*it);
     }
     cout << endl;
   } else {
     cout << "No Alpha sums have been specified " << endl;
   }
   
   if ( ioutfg.sumsb.size() != 0 ) {
     vector<int>::iterator it;
     cout << "Beta orbital sums are " << endl;
     for (it=(ioutfg.sumsb).begin();  it !=ioutfg.sumsb.end(); ++it){
       cout << " " << (*it);
     }
     cout << endl;
   } else {
     cout << "No Beta sums have been specified " << endl;
  }
}


/* A new method for printing out the set of determinants even if they 
   are disjoint distributed around nodes

   For now the format will probably result poor parallel output
   That is okay as this is really only used for debugging anyway
*/
void PsociDeterminants::printLocalSpatialsData()
{
#ifdef DETAILEDCHECK
  cout << GA::nodeid() << " Modified: printLocalSpatialsData" << endl;
  cout << GA::nodeid() << " Number of local dets is " << joutfg.size() << endl;
#endif

  int g_rank = GA::nodeid();
  int g_size = GA::nodes();
  for( int i=0; i< g_size; ++i ) {
    if ( i == g_rank ) {


  vector<JOUTFG>::iterator it;
  for(it=joutfg.begin(); it != joutfg.end(); ++it ) {
    JOUTFG ioutfg = (*it);
    cout << " completed  = " << ioutfg.completed_state << endl;
    cout << " index      = " << ioutfg.index << endl;
    cout << " nbf        = " << ioutfg.nbf << endl;
    cout << " jsym       = " << ioutfg.jsym << endl;
    cout << " ksym       = " << ioutfg.ksym << endl;
    cout << " mxopn      = " << ioutfg.mx_openshells << endl;
    cout << " nelec      = " << ioutfg.nelec << endl;
    //cout << " nalpha     = " << ioutfg.nalpha << endl;
    cout << " num_doubly = " << ioutfg.num_doubly << endl;
    cout << " num_singly = " << ioutfg.num_singly << endl;
    cout << " config: " << ioutfg.config << endl;
    cout << " nsefi      = " << ioutfg.nsefi << endl;
    cout << " ndeti      = " << ioutfg.ndeti << endl;
    vector<int>::iterator it;
    
    cout << "All reported orbital indexes have been shifted: A one was added. ";
    cout << " Internally, however, all orbitals begin at zero " << endl;
    
    cout << "List of doubly occupied orbs: " << endl;
    for ( it = ioutfg.doubly.begin(); it != ioutfg.doubly.end(); ++it ) {
      cout << " " << (*it)+1;
    }
    cout << endl << "List of singly occupied orbs: " << endl;
    for ( it = ioutfg.singly.begin(); it != ioutfg.singly.end(); ++it ) {
      cout << " " << (*it)+1;
    }
    cout << endl;
    
    if ( ioutfg.alpha_spins.size() != 0 ) {
      int ndet=0;
      vector<int>::iterator sit;
      sit = ioutfg.spin_projection.begin();
      cout << "List of alpha orbitals for this determinant " << endl;
      vector<vector<int> >::iterator vit;
      vector<int>::iterator inner;
      for (vit = ioutfg.alpha_spins.begin(); vit != ioutfg.alpha_spins.end(); ++vit) {
	++ndet;
	cout << "Determinant number: " << ndet <<" spin projection is " << (*sit) << " Number of alpha orbitals: " << (*vit).size()  <<  endl;
	for (inner=(*vit).begin(); inner !=(*vit).end(); ++inner){
	  cout << " " << (*inner)+1;
	}   
	cout << endl;
	++sit;
      } 
    } else {
      cout << "No alpha_spin orbitals have been specified " << endl;
    }
    
    if ( ioutfg.beta_spins.size() != 0 ) {
      int ndet=0;
      vector<int>::iterator sit;
      sit = ioutfg.spin_projection.begin();
      cout << "List of beta orbitals for this determinant " << endl;
      vector<vector<int> >::iterator vit;
      vector<int>::iterator inner;
      for (vit = ioutfg.beta_spins.begin(); vit != ioutfg.beta_spins.end(); ++vit) {
	++ndet;
	cout << "Determinant number: " << ndet <<" spin projection is " << (*sit)  << " Number of beta orbitals: " << (*vit).size()  <<  endl;
	for (inner=(*vit).begin(); inner !=(*vit).end(); ++inner){
	  cout << " " << (*inner)+1;
	} 
      cout << endl;
      ++sit;
      } 
    } else {
    cout << "No beta_spin orbitals have been specified " << endl;
    }
    
    if ( ioutfg.phases.size() != 0 ) {
      cout << "Determinant phases are " << endl;
      for (it=(ioutfg.phases).begin();  it !=ioutfg.phases.end(); ++it){
	cout << " " << (*it);
      }
      cout << endl;
    } else {
      cout << "No phases have been specified " << endl;
    }
    
    if ( ioutfg.sumsa.size() != 0 ) {
      //vector<int>::iterator it;
      cout << "Alpha orbital sums are " << endl;
      for (it=(ioutfg.sumsa).begin();  it !=ioutfg.sumsa.end(); ++it){
	cout << " " << (*it);
      }
      cout << endl;
    } else {
      cout << "No Alpha sums have been specified " << endl;
    }
    
    if ( ioutfg.sumsb.size() != 0 ) {
      vector<int>::iterator it;
      cout << "Beta orbital sums are " << endl;
      for (it=(ioutfg.sumsb).begin();  it !=ioutfg.sumsb.end(); ++it){
	cout << " " << (*it);
      }
      cout << endl;
    } else {
     cout << "No Beta sums have been specified " << endl;
    }
  }
    }
    GA::SERVICES.sync();
    //GA::sync();
  }

}
/* Generate all determinants with spin projection ( ms ) fractional needed to make
   double-group adapted functions from a spatial configuration. Formally named
   ( msfdet )
   
   was:
   call msfdet(iopn,jdbl,ms,ifigb,ifiga,-1)
*/

int PsociDeterminants::generateDeterminants( JOUTFG & ioutfg ) // msfdet
{
  int parity = fetchSpinParity(); //Applies to all dets
  int jsym = ioutfg.jsym; // spatial conf symmetry
  int ksym = l_ksym; // overall state symmetry
  short isgn = 1; // Triggers a spin flip (if Jsym even isgn = 1, else -1 )
  short iequal =-1; // Even/Odd options
  
  if ( (jsym + 1)/2 % 2 == 0 ) isgn = -1; //Jsym is ODD 
  if ( (jsym + 1)/2  == (ksym + 1)/2 ) iequal = 1; //Jsym = ksym
  
  if ( parity == 1 ) // E symmetry ( nel%2 = 1: odd -> fractional symmetry (1/2,3/2,etc)
    {
      if ( generateSpinFractionalDeterminants( ioutfg, isgn ) != 0 ) {
         cerr << "generateSpinFractionalDetermnants failed at " << endl;
         GA::Terminate();
      }
        } else {
      if ( iequal == 1 ) {
         if ( generateEVENSpinDeterminants( ioutfg ) != 0 ) {
         cerr << "generateEVENSpinDeterminants failed at " << endl;
         GA::Terminate();
         }
      } else {
         if ( generateODDSpinDeterminants( ioutfg ) != 0 ) {
         cerr << "generateODDSpinDeterminants failed at " << endl;
         GA::Terminate();
         }
      }
    }
#if 0
      cout << "size of alpha is " << ioutfg.alpha_spins.size() << endl;
      cout << "size of beta is " << ioutfg.beta_spins.size() << endl;
#endif
  
  return(0);
}

//Return the number of determinants generated ( should equal ndeti )
//Process ALL JOUTFG for the local set of configs

int PsociDeterminants::generateSpinFractionalDeterminants( JOUTFG & ioutfg, int isgn ) // msfdet
{
  if ( isgn != 1 && isgn != -1 ) {
    cerr << "Erroneous isgn was found: " << isgn << " Terminating " << endl;
    GA::Terminate();
  }
  
  int n_open = ioutfg.num_singly;
  
  vector<int> iopen = ioutfg.singly;
  vector<int> idbl = ioutfg.doubly;
  
#if 0
  printSpatialsData( ioutfg );
#endif
  
  int status = -1;
  switch ( n_open ) 
    {
    case PSOCI_ONE_OPEN:
      //cout << "One open shell electron found " << endl;
      status = oneOpenShellFractionalDeterminants( ioutfg, isgn );
      break;
    case PSOCI_THREE_OPEN:
      //cout << "Three open shell electrons found " << endl;
      status = threeOpenShellFractionalDeterminants( ioutfg, isgn );
      break;
    case PSOCI_FIVE_OPEN:
      //cout << "Five open shell electrons found " << endl;
      status = fiveOpenShellFractionalDeterminants( ioutfg, isgn );
      break;
    case PSOCI_SEVEN_OPEN:
      //cout << "Seven open shell electrons found " << endl;
      status = sevenOpenShellFractionalDeterminants( ioutfg, isgn );
      break;
    case PSOCI_NINE_OPEN:
      //cout << "nine open shell electrons found " << endl;
      status = nineOpenShellFractionalDeterminants( ioutfg, isgn );
      break;
    case PSOCI_ELEVEN_OPEN:
      //cout << "eleven open shell electrons found " << endl;
      status = elevenOpenShellFractionalDeterminants( ioutfg, isgn );
      break;
    default:
      cerr << "A maximum of  " << PSOCI_MAX_OPEN << " open shells is permitted: aborting " << endl;
      GA::Terminate();
    }
  
//Next for debugging only
  // printSpatialsData( ioutfg );
  return( status );
}

//Process ALL JOUTFG EVEN determinants for the local set of configs

int PsociDeterminants::generateEVENSpinDeterminants( JOUTFG & ioutfg ) // msedet
{
  
  int n_open = ioutfg.num_singly;
  
  vector<int> iopen = ioutfg.singly;
  vector<int> idbl = ioutfg.doubly;
  
#if 0
  printSpatialsData( ioutfg );
#endif
  
  int status = -1;
  switch ( n_open ) 
    {
    case PSOCI_ZERO_OPEN:
      //cout << "zero EVEN open shell electron found " << endl;
      status = zeroOpenShellEvenDeterminants( ioutfg );
      break;
    case PSOCI_TWO_OPEN:
      //cout << "two EVEN open shell electrons found " << endl;
      status = twoOpenShellEvenDeterminants( ioutfg );
      break;
    case PSOCI_FOUR_OPEN:
      //cout << "four EVEN open shell electrons found " << endl;
      status = fourOpenShellEvenDeterminants( ioutfg );
      break;
    case PSOCI_SIX_OPEN:
      //cout << "six EVEN open shell electrons found " << endl;
      status = sixOpenShellEvenDeterminants( ioutfg );
      break;
    case PSOCI_EIGHT_OPEN:
      //cout << "eight EVEN open shell electrons found " << endl;
      status = eightOpenShellEvenDeterminants( ioutfg );
      break;
    case PSOCI_TEN_OPEN:
      //cout << "ten EVEN open shell electrons found " << endl;
      status = tenOpenShellEvenDeterminants( ioutfg );
      break;
    default:
      cerr << "A maximum of  " << PSOCI_WHOLE_MAX_OPEN << " EVEN open shells is permitted: aborting " << endl;
      GA::Terminate();
    }
  // printSpatialsData( ioutfg );
  return( status );
}


//Process ALL JOUTFG ODD determinants for the local set of configs

int PsociDeterminants::generateODDSpinDeterminants( JOUTFG & ioutfg ) // msodet
{
  
  int n_open = ioutfg.num_singly;
  
  vector<int> iopen = ioutfg.singly;
  vector<int> idbl = ioutfg.doubly;
  
#if 0
  printSpatialsData( ioutfg );
#endif
  
  int status = -1;
  switch ( n_open ) 
    {
    case PSOCI_ZERO_OPEN:
      //cout << "zero ODD open shell electron found " << endl;
      /*
      cerr << "Zero electrons is an invalid call to generateODDSpinDeterminants:  at " << ioutfg.index << " Aborting " << endl;
      cerr << " config sym is " << ioutfg.jsym << " JOUTFG is " << ioutfg.config << endl;
      */
      status = zeroOpenShellOddDeterminants( ioutfg );
      break;
    case PSOCI_TWO_OPEN:
      //cout << "two ODD open shell electrons found " << endl;
      status = twoOpenShellOddDeterminants( ioutfg );
      break;
    case PSOCI_FOUR_OPEN:
      //cout << "four ODD open shell electrons found " << endl;
      status = fourOpenShellOddDeterminants( ioutfg );
      break;
    case PSOCI_SIX_OPEN:
      //cout << "six ODD open shell electrons found " << endl;
      status = sixOpenShellOddDeterminants( ioutfg );
      break;
    case PSOCI_EIGHT_OPEN:
      //cout << "eight ODD open shell electrons found " << endl;
      status = eightOpenShellOddDeterminants( ioutfg );
      break;
    case PSOCI_TEN_OPEN:
      //cout << "ten ODD open shell electrons found " << endl;
      status = tenOpenShellOddDeterminants( ioutfg );
      break;
    default:
      cerr << "A maximum of  " << PSOCI_WHOLE_MAX_OPEN << " ODD open shells is permitted: aborting " << endl;
      GA::Terminate();
    }
//   printSpatialsData( ioutfg );
  return( status );
}
/* Some variables are tersely named to maintain connection to the original code
   from which this is derived
*/

/* We use the same approach as that of Pitzer et al. Namely the alpha list contains ALL
   alphas ( including dbl occ alphas). Ditto for the beta. This seems wasteful in space
   but it will keep us close to their H construction approach and once inserted into the GA
   all *this space is terminated anyway.
*/

/* For the following sections all iteration bounds using Fortran indexing were retained. However, 
   For each variable, all indices were decremented by 1. This seemed like the most fool-proof way 
   to incorporate the existing code without adding a subtle bug...
*/

//The approach is cumbersome, but I really didn't want to significantly rewrite the determinant assembly code 

int PsociDeterminants::oneOpenShellFractionalDeterminants( JOUTFG & ioutfg, int isgn )
{
  int msbas = fetchMaxOpenShells() / 2; // For ALL spatial configurations
  int nalpha = ( l_nelec + 1 ) / 2; 
  int ndeti = ioutfg.ndeti;
  
  /*
    case 1:  one open-shell electron
    
    number of double-group functions = 1
    number of determinants           = 1
    
    the list of doubly occupied orbitals goes in both the spin-up
    and spin-down lists. the singly occupied orbital goes in the
    spin-up list.
    
  */
  
  if ( ioutfg.ndeti != 1 ) {
    cerr << "Inconsistent ndeti != 1 " << ioutfg.ndeti << " Aborting " << endl;
    GA::Terminate();
  }
  
  if ( ioutfg.singly.size() != 1 ) {
    cerr << "Inconsistent number of open shells" << ioutfg.singly.size();
    cerr << " configuration index " << ioutfg.index  << " Aborting " << endl;
    GA::Terminate();
  }
  
  vector<int> jdbl;
  jdbl.reserve( nalpha );
  
  vector< vector<int> > ifiga( ndeti, vector<int>(nalpha) );
  vector< vector<int> > ifigb( ndeti, vector<int>(nalpha) );
  vector<int> ms( ndeti );
  
  // Copy all current double occupied orbs
  vector<int>::iterator dit;
  
  for (dit = ioutfg.doubly.begin(); dit != ioutfg.doubly.end(); ++dit ) {
    jdbl.push_back( (*dit) );
  }
  
  int ndeta = 1;
  int num_doubly = ioutfg.doubly.size();
  ms[ndeta-1] = msbas + (1 + isgn)/2;
  
  ifigb[ndeta-1].resize(num_doubly);
  for( int jj=1; jj<= num_doubly; ++jj ) {
    ifigb[ndeta-1][jj-1] = jdbl[jj-1];
  }
  
  ifiga[ndeta-1].resize(nalpha);
  jdbl[nalpha-1] = ioutfg.singly[0];
  for( int jj=1; jj<= nalpha; ++jj ) {
    ifiga[ndeta-1][jj-1] = jdbl[jj-1];
  }
  
  if ( pushBackToJoutfg( ioutfg, ifiga, ifigb, ms, isgn ) != 0 ) {
    cerr << "Error: fatal error from pushBackToJoutfg: << Aborting " << endl;
    GA::Terminate();
  }
  return( 0 );
}   

// Next: 3 open shell orbitals
// Note: regarding the resizes...ndeta doesn't overlap the various code blocks
//      and only results < nalpha need resizing

int PsociDeterminants::threeOpenShellFractionalDeterminants( JOUTFG & ioutfg, int isgn )
{
  int msbas = fetchMaxOpenShells() / 2 ; // For ALL spatial configurations
  int nelec = l_nelec;
  int nalpha = ( nelec + 1 ) / 2; 
  
  /*
    case 2:  three open-shell electrons
    
    number of double-group functions = 4
    number of determinants           = 4
    
    arrangement of the singly occupied orbitals
    
    singly occupied orbitals
    det   alpha list   beta list   (sign) spin function
    
    1       1  2       3           (+) a a b
    2       1  3       2           (-) a b a
    3       2  3       1           (+) b a a
    4                  1  2  3     (+) b b b
    
  */
  if ( ioutfg.ndeti != 4 ) {
    cerr << "Inconsistent content of ndeti != 4 " << ioutfg.ndeti << " Aborting " << endl;
    GA::Terminate();
  }
  
  if ( ioutfg.singly.size() != 3 ) {
    cerr << "Inconsistent number of open shells" << ioutfg.singly.size();
    cerr << " configuration index " << ioutfg.index  << " Aborting " << endl;
    GA::Terminate();
  }
  
// Perform the bulk of the effort in local arrays
  
  vector<int> jdbl;
  
  int ndeti = ioutfg.ndeti;
  
  vector< vector<int> > ifiga( ndeti, vector<int>(nalpha) ); //  A reasonably starting size
  vector< vector<int> > ifigb( ndeti, vector<int>(nalpha) ); // ditto

  const int pad = 2;
  jdbl.reserve( nalpha + pad );

  vector<int> ms( ndeti ); 
  
  // Grab all current double occupied orbs
  vector<int>::iterator dit;
  for (dit = ioutfg.doubly.begin(); dit != ioutfg.doubly.end(); ++dit ) {
    jdbl.push_back( (*dit) );
  }
  int num_doubly = ioutfg.doubly.size();
  
  //   (alpha), then (beta), lists for determinants with ms = (1/2)
  
  int ndeta = 1;
  
  for(int i2 = 2; i2 <= 3; ++i2 ) {
    jdbl[nalpha-1] = ioutfg.singly[ i2-1 ];
    for( int i1 = 1 ; i1 <= i2-1; ++i1 ) {
      jdbl[nalpha-2] = ioutfg.singly[ i1-1 ];
      ms[ndeta-1] = (msbas + (1 + isgn)/2 );
      //ifiga[ndeta-1].resize(nalpha);
      for( int jj = 1; jj<= nalpha; ++jj ) {
	ifiga[ndeta-1][jj-1] = jdbl[jj-1];
      }
      ++ndeta;
    }
  }
  
  --nalpha;
  ndeta = 3; 
  for (int i1 = 1; i1 <= 3; ++i1 ) {
    jdbl[ nalpha-1] = ioutfg.singly[ i1-1 ];
    
    ifigb[ndeta-1].resize(nalpha);
    for (int jj=1; jj <= nalpha; ++jj ) {
      ifigb[ ndeta-1][jj-1] = jdbl[ jj-1 ];
    }
    --ndeta;
  }
  
  //   determinants with ms = (-3/2)
  
  ms[3] = msbas + (1 - 3*isgn)/2;
  ifiga[3].resize(num_doubly); 
  for( int jj=1; jj <=num_doubly; ++jj ) {
    ifiga[3][jj-1] = jdbl[jj-1];
  }
  
  for( int i1=1; i1<=3; ++i1) {
    jdbl[num_doubly+i1-1] = ioutfg.singly[ i1-1 ];
  }
  nalpha += 2;
  ifigb[3].resize(nalpha);
  for( int jj=1; jj <=nalpha; ++jj ) {
    ifigb[3][jj-1] = jdbl[jj-1];
  }
  
  //Push ms, ifiga, and ifigb onto the ioutfg vectors.
  
  if ( pushBackToJoutfg( ioutfg, ifiga, ifigb, ms, isgn ) != 0 ) {
    cerr << "Error: fatal error from pushBackToJoutfg: << Aborting " << endl;
    GA::Terminate();
  }
  return( 0 );
}

// Next 5 open shell orbitals
int PsociDeterminants::fiveOpenShellFractionalDeterminants( JOUTFG & ioutfg, int isgn )
{
  int msbas = fetchMaxOpenShells() / 2; // For ALL spatial configurations
  int nelec = l_nelec;
  int nalpha = ( nelec + 1 ) / 2; // Probably should make this a class variable instead of a struct member 
  
  /*
    case 3:  five open-shell electrons
    
    number of double-group functions = 16
    number of determinants           = 16
    
    arrangement of the singly occupied orbitals
    
    singly occupied orbitals
    det       alpha list        beta list     (sign) spin function
    
    1       1  2  3            4  5           (+) a a a b b
    2       1  2  4            3  5           (-) a a b a b
    3       1  3  4            2  5           (+) a b a a b
    4       2  3  4            1  5           (-) b a a a b
    5       1  2  5            3  4           (+) a a b b a
    6       1  3  5            2  4           (-) a b a b a
    7       2  3  5            1  4           (+) b a a b a
    8       1  4  5            2  3           (+) a b b a a
    9       2  4  5            1  3           (-) b a b a a
   10       3  4  5            1  2           (+) b b a a a
   11       5                  1  2  3  4     (+) b b b b a
   12       4                  1  2  3  5     (-) b b b a b
   13       3                  1  2  4  5     (+) b b a b b
   14       2                  1  3  4  5     (-) b a b b b
   15       1                  2  3  4  5     (+) a b b b b
   16       1  2  3  4  5                     (+) a a a a a
  */
  
 vector<int> jdbl;
 
 int ndeti = ioutfg.ndeti;
 
 vector< vector<int> > ifiga( ndeti, vector<int>(nalpha) );
 vector< vector<int> > ifigb( ndeti, vector<int>(nalpha) );

 const int pad = 2;
 jdbl.reserve( nalpha + pad );

 vector<int> ms( ndeti );
 
 vector<int>::iterator dit;
 for (dit = ioutfg.doubly.begin(); dit != ioutfg.doubly.end(); ++dit ) {
   jdbl.push_back( (*dit) );
 }
 
 int num_doubly = ioutfg.doubly.size();
 
 // (alpha), then (beta), lists for determinants with ms = (1/2)
 
 int ndeta = 1;
 for (int i3=3; i3 <= 5; ++i3 ) {
   jdbl[nalpha-1] = ioutfg.singly[i3-1]; 
   for (int i2=2; i2 <= i3-1; ++i2 ) {
     jdbl[nalpha-2] = ioutfg.singly[i2-1];
     for (int i1=1; i1 <= i2-1; ++i1) {
       jdbl[nalpha-3] = ioutfg.singly[i1-1];
       ms[ndeta-1] = (msbas + (1 + isgn)/2);
       
       //ifiga[ndeta-1].resize( nalpha );
       for (int jj=1; jj<=nalpha; ++jj) {
	 ifiga[ndeta-1][jj-1] = jdbl[jj-1];
       }
       ++ndeta;
     }
   }
 }
 
 --nalpha;
 ndeta = 10;
 for( int i2=2; i2 <=5; ++i2) {
   jdbl[nalpha-1] = ioutfg.singly[i2-1];
   for(int i1=1; i1 <= i2-1; ++i1) {
     jdbl[nalpha-2] = ioutfg.singly[i1-1];
     
     ifigb[ndeta-1].resize( nalpha );
     for (int jj=1; jj<=nalpha; ++jj) {
       ifigb[ndeta-1][jj-1] = jdbl[jj-1];
     }
     --ndeta;
   }
 }
 
 // (beta), then (alpha), lists for determinants with ms = (-3/2)
 
 nalpha += 2;
 ndeta = 11;
 for (int i4=4; i4<=5; ++i4) {
   jdbl[nalpha-1] = ioutfg.singly[i4-1];
   for (int i3=3; i3<=i4-1; ++i3) {
     jdbl[nalpha-2] = ioutfg.singly[i3-1];
     for (int i2=2; i2<=i3-1; ++i2) {
       jdbl[nalpha-3] = ioutfg.singly[i2-1];
       for (int i1=1; i1<=i2-1; ++i1) {
	 jdbl[nalpha-4] = ioutfg.singly[i1-1];
	 ms[ndeta-1] = (msbas + (1 - 3*isgn)/2);
         
	 ifigb[ndeta-1].resize( nalpha );
	 for(int jj=1; jj<= nalpha; ++jj) {
	   ifigb[ndeta-1][jj-1] = jdbl[jj-1];
	 }
	 ++ndeta;
       }
     }
   }
 }
 
 nalpha -= 3;
 ndeta = 15;
 for(int i1=1; i1<=5; ++i1) {
   jdbl[nalpha-1] = ioutfg.singly[i1-1];
   
   ifiga[ndeta-1].resize( nalpha );
   for(int jj=1; jj<= nalpha; ++jj) {
     ifiga[ndeta-1][jj-1] = jdbl[jj-1]; 
   }
   --ndeta;
 }
 
 // determinant with ms = (5/2)
 
    ms[15] = msbas + (1 + 5*isgn)/2;
    
    ifigb[15].resize( num_doubly ); // prune down empty zeros
    for(int jj=1; jj<=num_doubly; ++jj) {
      ifigb[15][jj-1] = jdbl[jj-1];
    }
    for(int i1=1; i1<=5; ++i1) {
      jdbl[num_doubly+i1-1] = ioutfg.singly[i1-1];
    }
    nalpha += 4;
    
    ifiga[15].resize( nalpha );
    for(int jj=1; jj<=nalpha; ++jj) {
      ifiga[15][jj-1] = jdbl[jj-1];
    }
    
    if ( pushBackToJoutfg( ioutfg, ifiga, ifigb, ms, isgn ) != 0 ) {
      cerr << "Error: fatal error from pushBackToJoutfg: << Aborting " << endl;
      GA::Terminate();
    }
    return( 0 );
} 

// Next 7 open shell orbitals
int PsociDeterminants::sevenOpenShellFractionalDeterminants( JOUTFG & ioutfg, int isgn )
{
  int msbas = fetchMaxOpenShells() / 2; // For ALL spatial configurations
  int nelec = l_nelec;
  int nalpha = ( nelec + 1 ) / 2; // Probably should make this a class variable instead of a struct member 
  
  /*
    case 4:  seven open-shell electrons
    
    number of double-group functions = 64 
    number of determinants           = 64 
    
  */
  
  vector<int> jdbl;
  
  int ndeti = ioutfg.ndeti;
  
  vector< vector<int> > ifiga( ndeti, vector<int>(nalpha) );
  vector< vector<int> > ifigb( ndeti, vector<int>(nalpha) );

  const int pad = 3;
  jdbl.reserve( nalpha + pad );

  vector<int> ms( ndeti );
  
  vector<int>::iterator dit;
  for (dit = ioutfg.doubly.begin(); dit != ioutfg.doubly.end(); ++dit ) {
    jdbl.push_back( (*dit) );
  }
  
  int num_doubly = ioutfg.doubly.size();
  
  //   (alpha), then (beta), lists for determinants with ms = (1/2)
  
  int ndeta = 1;
  for(int i4=4; i4<=7; ++i4) {
    jdbl[nalpha-1] = ioutfg.singly[i4-1];
    for(int i3=3; i3<=i4-1; ++i3) {
      jdbl[nalpha-2] = ioutfg.singly[i3-1];
      for(int i2=2; i2<=i3-1; ++i2) {
	jdbl[nalpha-3] = ioutfg.singly[i2-1];
	for(int i1=1; i1<=i2-1; ++i1) {
	  jdbl[nalpha-4] = ioutfg.singly[i1-1];
	  ms[ndeta-1] = (msbas + (1 + isgn)/2);
          
	  //ifiga[ndeta-1].resize( nalpha );
	  for (int jj=1; jj<=nalpha; ++jj) {
	    ifiga[ndeta-1][jj-1] = jdbl[jj-1];
	  }
	  ++ndeta;
	}
      }
    }
  }

  --nalpha;
  ndeta = 35;
  for(int i3=3;i3<=7; ++i3) {
    jdbl[nalpha-1] = ioutfg.singly[i3-1];
    for(int i2=2;i2<=i3-1; ++i2) {
      jdbl[nalpha-2] = ioutfg.singly[i2-1];
      for(int i1=1;i1<=i2-1; ++i1) { 
	jdbl[nalpha-3] = ioutfg.singly[i1-1];
	
	ifigb[ndeta-1].resize( nalpha );
	for (int jj=1; jj<=nalpha; ++jj) {
	  ifigb[ndeta-1][jj-1] = jdbl[jj-1];
	}
	--ndeta;
      }
    }
  }
  
  // (beta), then (alpha), lists for determinants with ms = (-3/2)
  
  nalpha += 2;
  ndeta = 36;
  for(int i5=5; i5<=7;++i5) {
    jdbl[nalpha-1] = ioutfg.singly[i5-1];
    for(int i4=4; i4<=i5-1;++i4) {
      jdbl[nalpha-2] = ioutfg.singly[i4-1];
      for(int i3=3; i3<=i4-1;++i3) {
	jdbl[nalpha-3] = ioutfg.singly[i3-1];
	for(int i2=2; i2<=i3-1;++i2) {
	  jdbl[nalpha-4] = ioutfg.singly[i2-1];
	  for(int i1=1; i1<=i2-1;++i1) {
	    jdbl[nalpha-5] = ioutfg.singly[i1-1];
	    ms[ndeta-1] = (msbas + (1 - 3*isgn)/2);
	    
	    ifigb[ndeta-1].resize( nalpha );
	    for (int jj=1; jj<=nalpha; ++jj) {
	      ifigb[ndeta-1][jj-1] = jdbl[jj-1];
	    }
	    ++ndeta;
	  }
	}
      }
    }
  }
  
  nalpha -= 3;
  ndeta = 56;
  
  for(int i2=2; i2<=7;++i2) {
    jdbl[nalpha-1] = ioutfg.singly[i2-1];
    for(int i1=1; i1<=i2-1;++i1) {
      jdbl[nalpha-2] = ioutfg.singly[i1-1];
      
      ifiga[ndeta-1].resize( nalpha );
      for (int jj=1; jj<=nalpha; ++jj) {
	ifiga[ndeta-1][jj-1] = jdbl[jj-1];
      }
      --ndeta;
    }
  }
  
  // (alpha), then (beta), lists for determinants with ms = (5/2)
  
  nalpha += 4;
  ndeta = 57;
  for(int i6=6; i6<=7; ++i6) {
    jdbl[nalpha-1] = ioutfg.singly[i6-1];
    for(int i5=5; i5<=i6-1;++i5) {
      jdbl[nalpha-2] = ioutfg.singly[i5-1];
      for(int i4=4; i4<=i5-1;++i4) {
	jdbl[nalpha-3] = ioutfg.singly[i4-1];
	for(int i3=3; i3<=i4-1;++i3) {
	  jdbl[nalpha-4] = ioutfg.singly[i3-1];
	  for(int i2=2; i2<=i3-1;++i2) {
	    jdbl[nalpha-5] = ioutfg.singly[i2-1];
	    for(int i1=1; i1<=i2-1;++i1) {
	      jdbl[nalpha-6] = ioutfg.singly[i1-1];
	      ms[ndeta-1] = (msbas + (1 + 5*isgn)/2);
	      
	      ifiga[ndeta-1].resize( nalpha );
	      for (int jj=1; jj<=nalpha; ++jj) {
		ifiga[ndeta-1][jj-1] = jdbl[jj-1];
	      }
	      ++ndeta;
	    }
	  }
	}
      }
    }
  }
  
  nalpha -= 5;
  ndeta = 63;
  for(int i1=1; i1<=7; ++i1) {
    jdbl[nalpha-1] = ioutfg.singly[i1-1];
    
    ifigb[ndeta-1].resize( nalpha );
    for (int jj=1; jj<=nalpha; ++jj) {
      ifigb[ndeta-1][jj-1] = jdbl[jj-1];
    }
    --ndeta;
  }
  
  // determinant with ms = (-7/2)
  
  ms[63] = msbas + (1 - 7*isgn)/2;
  ifiga[63].resize( num_doubly );
  for(int jj=1; jj<=num_doubly; ++jj) {
    ifiga[63][jj-1] = jdbl[jj-1];
  }
  for(int i1=1; i1<=7; ++i1) {
    jdbl[num_doubly+i1-1] = ioutfg.singly[i1-1];
  }
  nalpha += 6;
  
  ifigb[63].resize( nalpha );
  for(int jj=1; jj<=nalpha; ++jj) {
    ifigb[63][jj-1] = jdbl[jj-1];
  }
  
  if ( pushBackToJoutfg( ioutfg, ifiga, ifigb, ms, isgn ) != 0 ) {
    cerr << "Error: fatal error from pushBackToJoutfg: << Aborting " << endl;
    GA::Terminate();
  }   
  return( 0 );
}    

// Next 9 open shell orbitals
int PsociDeterminants::nineOpenShellFractionalDeterminants( JOUTFG & ioutfg, int isgn )
 {
   int msbas = fetchMaxOpenShells() / 2; // For ALL spatial configurations
   int nelec = l_nelec;
   int nalpha = ( nelec + 1 ) / 2; // Probably should make this a class variable instead of a struct member 
   
   /*
     case 5:  nine open-shell electrons
     
     number of double-group functions = 256 
     number of determinants           = 256 
     
   */
   
   vector<int> jdbl;
   
   int ndeti = ioutfg.ndeti;
   
   vector< vector<int> > ifiga( ndeti, vector<int>(nalpha) );
   vector< vector<int> > ifigb( ndeti, vector<int>(nalpha) );

   const int pad = 4;
   jdbl.reserve( nalpha + pad );

   vector<int> ms( ndeti );
   
   vector<int>::iterator dit;
   for (dit = ioutfg.doubly.begin(); dit != ioutfg.doubly.end(); ++dit ) {
     jdbl.push_back( (*dit) );
   }
   
   int num_doubly = ioutfg.doubly.size();
   
   //   (alpha), then (beta), lists for determinants with ms = (1/2)
   
   int ndeta = 1;
   for(int i5=5; i5<=9; ++i5) {
     jdbl[nalpha-1] = ioutfg.singly[i5-1];
     for(int i4=4; i4<=i5-1; ++i4) {
       jdbl[nalpha-2] = ioutfg.singly[i4-1];
       for(int i3=3; i3<=i4-1; ++i3) {
	 jdbl[nalpha-3] = ioutfg.singly[i3-1];
	 for(int i2=2; i2<=i3-1; ++i2) {
	   jdbl[nalpha-4] = ioutfg.singly[i2-1];
	   for(int i1=1; i1<=i2-1; ++i1) {
	     jdbl[nalpha-5] = ioutfg.singly[i1-1];
	     ms[ndeta-1] = (msbas + (1 + isgn)/2);
             
	     //ifiga[ndeta-1].resize( nalpha );
	     for (int jj=1; jj<=nalpha; ++jj) {
	       ifiga[ndeta-1][jj-1] = jdbl[jj-1];
	     }
             ++ndeta;
	   }
         }
       }
     }
   }
   
   --nalpha;
   ndeta = 126;
   for(int i4=4; i4<=9; ++i4) {
     jdbl[nalpha-1] = ioutfg.singly[i4-1];
     for(int i3=3;i3<=i4-1; ++i3) {
       jdbl[nalpha-2] = ioutfg.singly[i3-1];
       for(int i2=2;i2<=i3-1; ++i2) {
	 jdbl[nalpha-3] = ioutfg.singly[i2-1];
	 for(int i1=1;i1<=i2-1; ++i1) { 
           jdbl[nalpha-4] = ioutfg.singly[i1-1];
	   
           ifigb[ndeta-1].resize( nalpha );
	   for (int jj=1; jj<=nalpha; ++jj) {
	     ifigb[ndeta-1][jj-1] = jdbl[jj-1];
	   }
           --ndeta;
         }
       }
     }
   }
   
   // (beta), then (alpha), lists for determinants with ms = (-3/2)
   
   nalpha += 2;
   ndeta = 127;
   for(int i6=6; i6<=9;++i6) {
     jdbl[nalpha-1] = ioutfg.singly[i6-1];
     for(int i5=5; i5<=i6-1;++i5) {
       jdbl[nalpha-2] = ioutfg.singly[i5-1];
       for(int i4=4; i4<=i5-1;++i4) {
	 jdbl[nalpha-3] = ioutfg.singly[i4-1];
	 for(int i3=3; i3<=i4-1;++i3) {
	   jdbl[nalpha-4] = ioutfg.singly[i3-1];
	   for(int i2=2; i2<=i3-1;++i2) {
	     jdbl[nalpha-5] = ioutfg.singly[i2-1];
	     for(int i1=1; i1<=i2-1;++i1) {
	       jdbl[nalpha-6] = ioutfg.singly[i1-1];
	       ms[ndeta-1] = (msbas + (1 - 3*isgn)/2);
	       
	       ifigb[ndeta-1].resize( nalpha );
	       for (int jj=1; jj<=nalpha; ++jj) {
		 ifigb[ndeta-1][jj-1] = jdbl[jj-1];
	       }
	       ++ndeta;
	     }
	   }
	 }
       }
     }
   }
   
   nalpha -= 3;
   ndeta = 210;
   
   for(int i3=3; i3<=9; ++i3) {
     jdbl[nalpha-1] = ioutfg.singly[i3-1];
     for(int i2=2; i2<=i3-1;++i2) {
       jdbl[nalpha-2] = ioutfg.singly[i2-1];
       for(int i1=1; i1<=i2-1;++i1) {
	 jdbl[nalpha-3] = ioutfg.singly[i1-1];
	 
	 ifiga[ndeta-1].resize( nalpha );
	 for (int jj=1; jj<=nalpha; ++jj) {
	   ifiga[ndeta-1][jj-1] = jdbl[jj-1];
	 }
	 --ndeta;
       }
     }
   }
   
   // (alpha), then (beta), lists for determinants with ms = (5/2)
   
   nalpha += 4;
   ndeta = 211;
   for(int i7=7; i7<=9; ++i7) {
     jdbl[nalpha-1] = ioutfg.singly[i7-1];
     for(int i6=6; i6<=i7-1; ++i6) {
       jdbl[nalpha-2] = ioutfg.singly[i6-1];
       for(int i5=5; i5<=i6-1;++i5) {
	 jdbl[nalpha-3] = ioutfg.singly[i5-1];
	 for(int i4=4; i4<=i5-1;++i4) {
	   jdbl[nalpha-4] = ioutfg.singly[i4-1];
	   for(int i3=3; i3<=i4-1;++i3) {
             jdbl[nalpha-5] = ioutfg.singly[i3-1];
             for(int i2=2; i2<=i3-1;++i2) {
	       jdbl[nalpha-6] = ioutfg.singly[i2-1];
	       for(int i1=1; i1<=i2-1;++i1) {
		 jdbl[nalpha-7] = ioutfg.singly[i1-1];
		 ms[ndeta-1] = (msbas + (1 + 5*isgn)/2);
		 
		 ifiga[ndeta-1].resize( nalpha );
		 for (int jj=1; jj<=nalpha; ++jj) {
		   ifiga[ndeta-1][jj-1] = jdbl[jj-1];
		 }
		 ++ndeta;
	       }
	     }
           }
	 }
       }
     }
   }
   
   nalpha -= 5;
   ndeta = 246;
   for(int i2=2; i2<=9; ++i2) {
     jdbl[nalpha-1] = ioutfg.singly[i2-1];
     for(int i1=1; i1<=i2-1; ++i1) {
       jdbl[nalpha-2] = ioutfg.singly[i1-1];
       
       ifigb[ndeta-1].resize( nalpha );
       for (int jj=1; jj<=nalpha; ++jj) {
	 ifigb[ndeta-1][jj-1] = jdbl[jj-1];
       }
       --ndeta;
     }
   }
   
   // determinant with ms = (-7/2)
 
   nalpha += 6;
   ndeta = 247;
   for(int i8=8; i8<=9; ++i8) {
     jdbl[nalpha-1] = ioutfg.singly[i8-1];
     for(int i7=7; i7<=i8-1; ++i7) {
       jdbl[nalpha-2] = ioutfg.singly[i7-1];
       for(int i6=6; i6<=i7-1; ++i6) {
	 jdbl[nalpha-3] = ioutfg.singly[i6-1];
	 for(int i5=5; i5<=i6-1;++i5) {
	   jdbl[nalpha-4] = ioutfg.singly[i5-1];
	   for(int i4=4; i4<=i5-1;++i4) {
	     jdbl[nalpha-5] = ioutfg.singly[i4-1];
	     for(int i3=3; i3<=i4-1;++i3) {
	       jdbl[nalpha-6] = ioutfg.singly[i3-1];
	       for(int i2=2; i2<=i3-1;++i2) {
		 jdbl[nalpha-7] = ioutfg.singly[i2-1];
		 for(int i1=1; i1<=i2-1;++i1) {
		   jdbl[nalpha-8] = ioutfg.singly[i1-1];
		   ms[ndeta-1] = (msbas + (1 - 7*isgn)/2);
		   
		   ifigb[ndeta-1].resize( nalpha );
		   for (int jj=1; jj<=nalpha; ++jj) {
		     ifigb[ndeta-1][jj-1] = jdbl[jj-1];
		   }
                   ++ndeta;
                   }
               }
	     }
	   }
	 }
       }
     }
   }
   
   nalpha -= 7;
   ndeta = 255;
   for(int i1=1;i1<=9;++i1) {
     jdbl[nalpha-1] = ioutfg.singly[i1-1];
     
     ifiga[ndeta-1].resize( nalpha );
     for (int jj=1; jj<=nalpha; ++jj) {
       ifiga[ndeta-1][jj-1] = jdbl[jj-1];
     }
     --ndeta;
   }
   
   // determinant with ms = (9/2)
   
   ms[255] = msbas + (1 + 9*isgn)/2;
   ifigb[255].resize( num_doubly );
   for(int jj=1; jj<=num_doubly; ++jj) {
     ifigb[255][jj-1] = jdbl[jj-1];
   }
   for(int i1=1; i1<=9; ++i1) {
     jdbl[num_doubly+i1-1] = ioutfg.singly[i1-1];
   }
   nalpha += 8;
   
   ifiga[255].resize( nalpha );
   for(int jj=1; jj<=nalpha; ++jj) {
     ifiga[255][jj-1] = jdbl[jj-1];
   }
   
   if ( pushBackToJoutfg( ioutfg, ifiga, ifigb, ms, isgn ) != 0 ) {
     cerr << "Error: fatal error from pushBackToJoutfg: << Aborting " << endl;
     GA::Terminate();
   }   
   return( 0 );
 }    

// Next 11 open shell orbitals
int PsociDeterminants::elevenOpenShellFractionalDeterminants( JOUTFG & ioutfg, int isgn )
{
  int msbas = fetchMaxOpenShells() / 2; // For ALL spatial configurations
  int nelec = l_nelec;
  int nalpha = ( nelec + 1 ) / 2; // Probably should make this a class variable instead of a struct member 
  
  /*
    case 5:  eleven open-shell electrons
    
    number of double-group functions = 1024
    number of determinants           = 1024
    
  */
  
  vector<int> jdbl;
  
  int ndeti = ioutfg.ndeti;
  
  const int pad = 5;
  
  vector< vector<int> > ifiga( ndeti, vector<int>(nalpha) );
  vector< vector<int> > ifigb( ndeti, vector<int>(nalpha) );
  jdbl.reserve( nalpha + pad );
  
  vector<int> ms( ndeti );
  
  vector<int>::iterator dit;
  for (dit = ioutfg.doubly.begin(); dit != ioutfg.doubly.end(); ++dit ) {
    jdbl.push_back( (*dit) );
  }
  
  int num_doubly = ioutfg.doubly.size();
  
  //   (alpha), then (beta), lists for determinants with ms = (1/2)
  
  int ndeta = 1;
  for(int i6=6; i6<=11; ++i6) {
    jdbl[nalpha-1] = ioutfg.singly[i6-1];
    for(int i5=5; i5<=i6-1; ++i5) {
      jdbl[nalpha-2] = ioutfg.singly[i5-1];
      for(int i4=4; i4<=i5-1; ++i4) {
	jdbl[nalpha-3] = ioutfg.singly[i4-1];
	for(int i3=3; i3<=i4-1; ++i3) {
	  jdbl[nalpha-4] = ioutfg.singly[i3-1];
	  for(int i2=2; i2<=i3-1; ++i2) {
	    jdbl[nalpha-5] = ioutfg.singly[i2-1];
	    for(int i1=1; i1<=i2-1; ++i1) {
	      jdbl[nalpha-6] = ioutfg.singly[i1-1];
	      ms[ndeta-1] = (msbas + (1 + isgn)/2);
	      
	      //ifiga[ndeta-1].resize( nalpha );
	      for (int jj=1; jj<=nalpha; ++jj) {
		ifiga[ndeta-1][jj-1] = jdbl[jj-1];
	      }
	      ++ndeta;
	    }
	  }
	}
      }
    }
  }
  
  --nalpha;
  ndeta = 462;
  for(int i5=5; i5<=11; ++i5) {
    jdbl[nalpha-1] = ioutfg.singly[i5-1];
    for(int i4=4; i4<=i5-1; ++i4) {
      jdbl[nalpha-2] = ioutfg.singly[i4-1];
      for(int i3=3;i3<=i4-1; ++i3) {
	jdbl[nalpha-3] = ioutfg.singly[i3-1];
	for(int i2=2;i2<=i3-1; ++i2) {
	  jdbl[nalpha-4] = ioutfg.singly[i2-1];
	  for(int i1=1;i1<=i2-1; ++i1) { 
	    jdbl[nalpha-5] = ioutfg.singly[i1-1];
	    
	    ifigb[ndeta-1].resize( nalpha );
	    for (int jj=1; jj<=nalpha; ++jj) {
	      ifigb[ndeta-1][jj-1] = jdbl[jj-1];
	    }
	    --ndeta;
	  }
	}
      }
    }
  }
  
  // (beta), then (alpha), lists for determinants with ms = (-3/2)
  
  nalpha += 2;
  ndeta = 463;
  for(int i7=7; i7<=11;++i7) {
    jdbl[nalpha-1] = ioutfg.singly[i7-1];
    for(int i6=6; i6<=i7-1;++i6) {
      jdbl[nalpha-2] = ioutfg.singly[i6-1];
      for(int i5=5; i5<=i6-1;++i5) {
	jdbl[nalpha-3] = ioutfg.singly[i5-1];
	for(int i4=4; i4<=i5-1;++i4) {
	  jdbl[nalpha-4] = ioutfg.singly[i4-1];
	  for(int i3=3; i3<=i4-1;++i3) {
	    jdbl[nalpha-5] = ioutfg.singly[i3-1];
	    for(int i2=2; i2<=i3-1;++i2) {
	      jdbl[nalpha-6] = ioutfg.singly[i2-1];
	      for(int i1=1; i1<=i2-1;++i1) {
		jdbl[nalpha-7] = ioutfg.singly[i1-1];
		ms[ndeta-1] = (msbas + (1 - 3*isgn)/2);
		
		ifigb[ndeta-1].resize( nalpha );
		for (int jj=1; jj<=nalpha; ++jj) {
		  ifigb[ndeta-1][jj-1] = jdbl[jj-1];
		}
		++ndeta;
	      }
	    }
	  }
	}
      }
    }
  }
  
  nalpha -= 3;
  ndeta = 792;
  
  for(int i4=4; i4<=11; ++i4) {
    jdbl[nalpha-1] = ioutfg.singly[i4-1];
    
    for(int i3=3; i3<=i4-1; ++i3) {
      jdbl[nalpha-2] = ioutfg.singly[i3-1];
      for(int i2=2; i2<=i3-1;++i2) {
	jdbl[nalpha-3] = ioutfg.singly[i2-1];
	for(int i1=1; i1<=i2-1;++i1) {
	  jdbl[nalpha-4] = ioutfg.singly[i1-1];
	  
	  ifiga[ndeta-1].resize( nalpha );
	  for (int jj=1; jj<=nalpha; ++jj) {
	    ifiga[ndeta-1][jj-1] = jdbl[jj-1];
	  }
	  --ndeta;
	}
      }
    }
  }
  
  // (alpha), then (beta), lists for determinants with ms = (5/2)
   
  nalpha += 4;
  ndeta = 793;
  for(int i8=8; i8<=11; ++i8) {
    jdbl[nalpha-1] = ioutfg.singly[i8-1];
    
    for(int i7=7; i7<=i8-1; ++i7) {
      jdbl[nalpha-2] = ioutfg.singly[i7-1];
      for(int i6=6; i6<=i7-1; ++i6) {
	jdbl[nalpha-3] = ioutfg.singly[i6-1];
	for(int i5=5; i5<=i6-1;++i5) {
	  jdbl[nalpha-4] = ioutfg.singly[i5-1];
	  for(int i4=4; i4<=i5-1;++i4) {
	    jdbl[nalpha-5] = ioutfg.singly[i4-1];
	    for(int i3=3; i3<=i4-1;++i3) {
	      jdbl[nalpha-6] = ioutfg.singly[i3-1];
	      for(int i2=2; i2<=i3-1;++i2) {
		jdbl[nalpha-7] = ioutfg.singly[i2-1];
		for(int i1=1; i1<=i2-1;++i1) {
		  jdbl[nalpha-8] = ioutfg.singly[i1-1];
		  ms[ndeta-1] = (msbas + (1 + 5*isgn)/2);
		  
		  ifiga[ndeta-1].resize( nalpha );
		  for (int jj=1; jj<=nalpha; ++jj) {
		    ifiga[ndeta-1][jj-1] = jdbl[jj-1];
		  }
		  ++ndeta;
		}
	      }
	    }
	  }
	}
      }
    }
  }
  
  nalpha -= 5;
  ndeta = 957;
  for(int i3=3; i3<=11; ++i3) {
    jdbl[nalpha-1] = ioutfg.singly[i3-1];
    
    for(int i2=2; i2<=i3-1; ++i2) {
      jdbl[nalpha-2] = ioutfg.singly[i2-1];
      for(int i1=1; i1<=i2-1; ++i1) {
	jdbl[nalpha-3] = ioutfg.singly[i1-1];
	
	ifigb[ndeta-1].resize( nalpha );
	for (int jj=1; jj<=nalpha; ++jj) {
	  ifigb[ndeta-1][jj-1] = jdbl[jj-1];
	}
	--ndeta;
      }
    }
  }
  
  // determinant with ms = (-7/2)
  
   nalpha += 6;
   ndeta = 958;

   for(int i9=9; i9<=11; ++i9) {
     jdbl[nalpha-1] = ioutfg.singly[i9-1];
     for(int i8=8; i8<=i9-1; ++i8) {
       jdbl[nalpha-2] = ioutfg.singly[i8-1];
       for(int i7=7; i7<=i8-1; ++i7) {
	 jdbl[nalpha-3] = ioutfg.singly[i7-1];
	 for(int i6=6; i6<=i7-1; ++i6) {
	   jdbl[nalpha-4] = ioutfg.singly[i6-1];
	   for(int i5=5; i5<=i6-1;++i5) {
	     jdbl[nalpha-5] = ioutfg.singly[i5-1];
	     for(int i4=4; i4<=i5-1;++i4) {
	       jdbl[nalpha-6] = ioutfg.singly[i4-1];
	       for(int i3=3; i3<=i4-1;++i3) {
		 jdbl[nalpha-7] = ioutfg.singly[i3-1];
		 for(int i2=2; i2<=i3-1;++i2) {
		   jdbl[nalpha-8] = ioutfg.singly[i2-1];
		   for(int i1=1; i1<=i2-1;++i1) {
		     jdbl[nalpha-9] = ioutfg.singly[i1-1];
		     ms[ndeta-1] = (msbas + (1 - 7*isgn)/2);
		     
		     ifigb[ndeta-1].resize( nalpha );
		     for (int jj=1; jj<=nalpha; ++jj) {
		       ifigb[ndeta-1][jj-1] = jdbl[jj-1];
		     }
		     ++ndeta;
                   }
		 }
	       }
	     }
	   }
	 }
       }
     }
   }
   
   nalpha -= 7;
   ndeta = 1012;
   for(int i2=2;i2<=11;++i2) {
     jdbl[nalpha-1] = ioutfg.singly[i2-1];
     for(int i1=1;i1<=i2-1;++i1) {
       jdbl[nalpha-2] = ioutfg.singly[i1-1];

       ifiga[ndeta-1].resize( nalpha );
       for (int jj=1; jj<=nalpha; ++jj) {
	 ifiga[ndeta-1][jj-1] = jdbl[jj-1];
       }
       --ndeta;
     }
   }
   
   // (alpha), then (beta), lists for determinants with ms = (9/2)

   nalpha += 8;
   ndeta = 1013;
   
   for(int i10=10; i10<=11; ++i10) {
     jdbl[nalpha-1] = ioutfg.singly[i10-1];
     for(int i9=9; i9<=i10-1; ++i9) {
       jdbl[nalpha-2] = ioutfg.singly[i9-1];
       for(int i8=8; i8<=i9-1; ++i8) {
	 jdbl[nalpha-3] = ioutfg.singly[i8-1];
	 for(int i7=7; i7<=i8-1; ++i7) {
	   jdbl[nalpha-4] = ioutfg.singly[i7-1];
	   for(int i6=6; i6<=i7-1; ++i6) {
	     jdbl[nalpha-5] = ioutfg.singly[i6-1];
	     for(int i5=5; i5<=i6-1;++i5) {
	       jdbl[nalpha-6] = ioutfg.singly[i5-1];
	       for(int i4=4; i4<=i5-1;++i4) {
		 jdbl[nalpha-7] = ioutfg.singly[i4-1];
		 for(int i3=3; i3<=i4-1;++i3) {
		   jdbl[nalpha-8] = ioutfg.singly[i3-1];
		   for(int i2=2; i2<=i3-1;++i2) {
		     jdbl[nalpha-9] = ioutfg.singly[i2-1];
		     for(int i1=1; i1<=i2-1;++i1) {
		       jdbl[nalpha-10] = ioutfg.singly[i1-1];
		       ms[ndeta-1] = (msbas + (1 + 9*isgn)/2);
		       
		       ifiga[ndeta-1].resize( nalpha );
		       for (int jj=1; jj<=nalpha; ++jj) {
			 ifiga[ndeta-1][jj-1] = jdbl[jj-1];
		       }
		       ++ndeta;
		     }
		   }
		 }
	       }
	     }
	   }
	 }
       }
     }
   }
   
   nalpha -= 9;
   ndeta = 1023;
   
   for(int i1=1;i1<=11;++i1) {
     jdbl[nalpha-1] = ioutfg.singly[i1-1];
     
     ifigb[ndeta-1].resize( nalpha );
     for (int jj=1; jj<=nalpha; ++jj) {
       ifigb[ndeta-1][jj-1] = jdbl[jj-1];
     }
     --ndeta;
   }
   
   // determinant with ms = (-11/2)
   
   ms[1023] = msbas + (1 - 11*isgn)/2;
   ifiga[1023].resize( num_doubly );
   for(int jj=1; jj<=num_doubly; ++jj) {
     ifiga[1023][jj-1] = jdbl[jj-1];
   }
   for(int i1=1; i1<=11; ++i1) {
     jdbl[num_doubly+i1-1] = ioutfg.singly[i1-1];
   }
   nalpha += 10;
   
   ifigb[1023].resize( nalpha );
   for(int jj=1; jj<=nalpha; ++jj) {
     ifigb[1023][jj-1] = jdbl[jj-1];
   }
   
   if ( pushBackToJoutfg( ioutfg, ifiga, ifigb, ms, isgn ) != 0 ) {
     cerr << "Error: fatal error from pushBackToJoutfg: << Aborting " << endl;
     GA::Terminate();
   }   
   return( 0 );
}    

// Local: Operates only on the local Determinants data  - Strong implied formatting issues - conforms with original code
// Not recommended to make this public
int PsociDeterminants::pushBackToJoutfg( JOUTFG & ioutfg, vector< vector<int> > & ifiga , vector< vector<int> > & ifigb, vector<int> &  ms )
{
  
  int ndeti = ioutfg.ndeti;
  if ( ifiga.size () != ndeti || ifigb.size () != ndeti || ms.size() != ndeti ) {
    cerr << "Unexpected dimensions for ifiga, ifigb, or ms: size is too big: " << endl;
    cerr << " ifiga = " << ifiga.size () << " ifigb = " << ifigb.size () << " ms = " << ms.size() << endl;
    GA::Terminate();
  }
  
  // Arrays behave as A[ndeti][nalpha] == vector<vector<int> > A(ndeti, vector<int>(nalpha) )
  // The column-major (implied) order is retained  to minimize changes to original code. Most
  // of the individual data objects fit in L1/L2 anyway.
  
  
  ioutfg.alpha_spins = ifiga;
  ioutfg.beta_spins = ifigb;
  ioutfg.spin_projection = ms;
  
/* New stuff Dec 2, 2011 Construct values for nalpha and nbeta for each orbital list
   Then populate a conforming array and add to the new joutfg members. This is required
   because in constructing H, we do not want to be constantly resizing alpha/beta arrays 
*/

/* This will work because for the determinant creation steps all arrays are of correct length
   This is not true for the H construction, however, where the temp arrays are of MAXIMUM length
   to save time with .resize() operations

   As for storage to GA, we already stored them
*/

   ioutfg.nalpha.clear();
   ioutfg.nbeta.clear();

   vector<vector<int> >::iterator it;
   for( it =ifiga.begin(); it != ifiga.end(); ++it ) {
      ioutfg.nalpha.push_back( (*it).size() );
   }

   for( it =ifigb.begin(); it != ifigb.end(); ++it ) {
      ioutfg.nbeta.push_back( (*it).size() );
   }

#ifdef DETAILEDCHECK
  cout << "size of new nalpha is " << ioutfg.nalpha.size();
  cout << "size of new nbeta is " << ioutfg.nbeta.size();
#endif

  return(0);
}

/* Used only for fractional determinantes since depending on isign we must swap alphas and betas
*/
int PsociDeterminants::pushBackToJoutfg( JOUTFG & ioutfg, vector< vector<int> > & ifiga , vector< vector<int> > & ifigb, vector<int> &  ms, int isign )
{
    if ( isign == -1 ) {
       return( pushBackToJoutfg( ioutfg, ifigb , ifiga, ms) );
    } else if ( isign == 1 ) {
       return( pushBackToJoutfg( ioutfg, ifiga , ifigb,  ms) );
    } else {
       GA::error(" wrong isign value ", isign );
    }
 return(-1);
}


// Begin what formally was MSEDET generators- even symmetry determinants

int PsociDeterminants::zeroOpenShellEvenDeterminants( JOUTFG & ioutfg )
{
  int msbas = fetchMaxOpenShells() / 2; 
  // int nalpha =  l_nelec / 2; 
  int ndeti = ioutfg.ndeti;
  
  /*
    case 1:  zero open-shell electron
    
    number of double-group functions = 1
    number of determinants           = 1
    
    the list of doubly occupied orbitals goes in both the spin-up
    and spin-down lists. the singly occupied orbital goes in the
    spin-up list.
    
  */
  
  if ( ioutfg.ndeti != 1 ) {
    cerr << "Inconsistent ndeti != 1 " << ioutfg.ndeti << " Aborting " << endl;
    GA::Terminate();
  }
  
  if ( ioutfg.singly.size() != 0 ) {
    cerr << "Inconsistent number of open shells" << ioutfg.singly.size();
    cerr << " configuration index " << ioutfg.index  << " Aborting " << endl;
    GA::Terminate();
  }
  
  int num_doubly = ioutfg.doubly.size();
  
  vector<int> jdbl;
  jdbl.reserve( num_doubly );
  
  vector< vector<int> > ifiga( ndeti, vector<int>(num_doubly) );
  vector< vector<int> > ifigb( ndeti, vector<int>(num_doubly) );
  vector<int> ms( ndeti );
  
  // Copy all current double occupied orbs
  vector<int>::iterator dit;
  
  for (dit = ioutfg.doubly.begin(); dit != ioutfg.doubly.end(); ++dit ) {
    jdbl.push_back( (*dit) );
  }
  
  int ndeta = 1;
  ms[0] = msbas;
  
  ifiga[ndeta-1].resize(num_doubly);
  ifigb[ndeta-1].resize(num_doubly);
  for( int jj=1; jj<= num_doubly; ++jj ) {
    ifiga[ndeta-1][jj-1] = jdbl[jj-1];
    ifigb[ndeta-1][jj-1] = jdbl[jj-1];
  }
  
  if ( pushBackToJoutfg( ioutfg, ifiga, ifigb, ms ) != 0 ) {
    cerr << "Error: fatal error from pushBackToJoutfg: << Aborting " << endl;
    GA::Terminate();
  }
  return( 0 );
}   

int PsociDeterminants::twoOpenShellEvenDeterminants( JOUTFG & ioutfg )
{
  int msbas = fetchMaxOpenShells() / 2; 
  int nalpha = l_nelec / 2; 
  int ndeti = ioutfg.ndeti;
  int ndeta = 1;
  
  /*
    case 2:  two open-shell electrons
    
    number of double-group functions = 1
    number of determinants           = 2
    
    arrangement of the singly occupied orbitals
    
    singly occupied orbitals
    det  alpha list  beta list   (sign) spin function
    
    1        1          2        (+) a b
    2        2          1        (-) b a
    
    the spin function and its sign are as obtained when the singly
    occupied orbitals are arranged in ascending order
    
    the doubly occupied orbitals are given first in the lists
    
  */
  
  if ( ioutfg.ndeti != 2 ) {
    cerr << "Inconsistent ndeti != 2 " << ioutfg.ndeti << " Aborting " << endl;
    GA::Terminate();
  }
  
  if ( ioutfg.singly.size() != 2 ) {
    cerr << "Inconsistent number of open shells" << ioutfg.singly.size();
    cerr << " configuration index " << ioutfg.index  << " Aborting " << endl;
    GA::Terminate();
  }
  
  vector<int> jdbl;
  jdbl.reserve( nalpha );
  
  vector< vector<int> > ifiga( ndeti, vector<int>(nalpha) );
  vector< vector<int> > ifigb( ndeti, vector<int>(nalpha) );
  vector<int> ms( ndeti );
  
  vector<int>::iterator dit;
  
  for (dit = ioutfg.doubly.begin(); dit != ioutfg.doubly.end(); ++dit ) {
    jdbl.push_back( (*dit) );
  }
  
  ms[0] = msbas; 
  ms[1] = msbas;
  
  
  jdbl[nalpha-1] = ioutfg.singly[0];
  for( int jj=1; jj<= nalpha; ++jj ) {
    ifiga[ndeta-1][jj-1] = jdbl[jj-1];
    ifigb[ndeta][jj-1] = jdbl[jj-1];
  }
  
  //ifiga[ndeta].resize(nalpha);
  //ifigb[ndeta-1].resize(nalpha);
  jdbl[nalpha-1] = ioutfg.singly[1];
  for( int jj=1; jj<= nalpha; ++jj ) {
    ifiga[ndeta][jj-1] = jdbl[jj-1];
    ifigb[ndeta-1][jj-1] = jdbl[jj-1];
  }
  
  if ( pushBackToJoutfg( ioutfg, ifiga, ifigb, ms ) != 0 ) {
    cerr << "Error: fatal error from pushBackToJoutfg: << Aborting " << endl;
    GA::Terminate();
  }
  return( 0 );
}   

//Four electron EVEN  problem
int PsociDeterminants::fourOpenShellEvenDeterminants( JOUTFG & ioutfg )
{
  int msbas = fetchMaxOpenShells() / 2;
  int nalpha = l_nelec / 2;
  int ndeti = ioutfg.ndeti;
  
  /*
    case 3:  four open-shell electrons
    
    number of double-group functions = 4
    number of determinants           = 8
    
    arrangement of the singly occupied orbitals
    
    singly occupied orbitals
    det      alpha list    beta list    (sign) spin function
    
    1       1  2          3  4          (+) a a b b
    2       3  4          1  2          (+) b b a a
    3       1  3          2  4          (-) a b a b
    4       2  4          1  3          (-) b a b a
    5       2  3          1  4          (+) b a a b
    6       1  4          2  3          (+) a b b a
    7       1  2  3  4                  (+) a a a a
    8                     1  2  3  4    (+) b b b b
  */
  
  if ( ioutfg.ndeti != 8 ) {
    cerr << "Inconsistent ndeti != 1 " << ioutfg.ndeti << " Aborting " << endl;
    GA::Terminate();
  } 
  
  if ( ioutfg.singly.size() != 4 ) {
    cerr << "Inconsistent number of open shells" << ioutfg.singly.size();
    cerr << " configuration index " << ioutfg.index  << " Aborting " << endl;
    GA::Terminate();
  } 
  
  const int pad = 2;
  vector<int> jdbl;
  jdbl.reserve( nalpha + pad );
  
  vector< vector<int> > ifiga( ndeti, vector<int>(nalpha) );
  vector< vector<int> > ifigb( ndeti, vector<int>(nalpha) );
  vector<int> ms( ndeti );
  
  vector<int>::iterator dit; 
  
  for (dit = ioutfg.doubly.begin(); dit != ioutfg.doubly.end(); ++dit ) {
    jdbl.push_back( (*dit) );
  }
  
  // determinants with ms = 0
  
  int ndeta = 1;
  int ndetb, ndetc;
  int num_doubly = ioutfg.doubly.size();
  
  for( int i2=2; i2<=4; ++i2) {
    jdbl[nalpha-1] = ioutfg.singly[i2-1];
    for( int i1=1; i1<=i2-1; ++i1) {
      jdbl[nalpha-2] = ioutfg.singly[i1-1];
      
      ndetb = min( ndeta, 13-ndeta );
      ndetc = min( ndeta+1, 12-ndeta);
      ms[ndetb-1] = msbas;
      
      //ifiga[ndetb-1].resize(nalpha);
      //ifigb[ndetc-1].resize(nalpha);
      for( int jj=1; jj<= nalpha; ++jj ) {
	ifiga[ndetb-1][jj-1] = jdbl[jj-1];
	ifigb[ndetc-1][jj-1] = jdbl[jj-1];
      }
      ndeta += 2;
    }
  }
  
  // determinants with ms = 2, -2
  
  ms[6] = msbas + 2;
  ms[7] = msbas - 2;
  
  ifiga[7].resize(num_doubly);
  ifigb[6].resize(num_doubly);
  for( int jj=1; jj<= num_doubly; ++jj ) {
    ifiga[7][jj-1] = jdbl[jj-1];
    ifigb[6][jj-1] = jdbl[jj-1];
  }
  for( int i1=1; i1<=4; ++i1) {
    jdbl[num_doubly+i1-1] = ioutfg.singly[i1-1];
  }
  
  nalpha += 2;
  ifiga[6].resize(nalpha);
  ifigb[7].resize(nalpha);
  for( int jj=1; jj<= nalpha; ++jj ) {
    ifiga[6][jj-1] = jdbl[jj-1];
    ifigb[7][jj-1] = jdbl[jj-1];
  }
  
  if ( pushBackToJoutfg( ioutfg, ifiga, ifigb, ms ) != 0 ) {
    cerr << "Error: fatal error from pushBackToJoutfg: << Aborting " << endl;
    GA::Terminate();
  }
  return( 0 );
}

// six electrons EVEN problem
int PsociDeterminants::sixOpenShellEvenDeterminants( JOUTFG & ioutfg )
{
  int msbas = fetchMaxOpenShells() / 2;
  int nalpha = l_nelec / 2;
  int ndeti = ioutfg.ndeti;
  
 /*
  number of double-group functions = 16
  number of determinants           = 32
  
  arrangement of the singly occupied orbitals
  
  singly occupied orbitals
  det       alpha list        beta list     (sign) spin function
  
   1       1  2  3          4  5  6          (+) a a a b b b
   2       4  5  6          1  2  3          (-) b b b a a a
   3       1  2  4          3  5  6          (-) a a b a b b
   4       3  5  6          1  2  4          (+) b b a b a a
   5       1  3  4          2  5  6          (+) a b a a b b
   6       2  5  6          1  3  4          (-) b a b b a a
   7       2  3  4          1  5  6          (-) b a a a b b
   8       1  5  6          2  3  4          (+) a b b b a a
   9       1  2  5          3  4  6          (+) a a b b a b
  10       3  4  6          1  2  5          (-) b b a a b a
  11       1  3  5          2  4  6          (-) a b a b a b
  12       2  4  6          1  3  5          (+) b a b a b a
  13       2  3  5          1  4  6          (+) b a a b a b
  14       1  4  6          2  3  5          (-) a b b a b a
  15       1  4  5          2  3  6          (+) a b b a a b
  16       2  3  6          1  4  5          (-) b a a b b a
  17       2  4  5          1  3  6          (-) b a b a a b
  18       1  3  6          2  4  5          (+) a b a b b a
  19       3  4  5          1  2  6          (+) b b a a a b
  20       1  2  6          3  4  5          (-) a a b b b a
  21       1  2  3  4  5    6                (+) a a a a a b
  22       6                1  2  3  4  5    (-) b b b b b a
  23       1  2  3  4  6    5                (-) a a a a b a
  24       5                1  2  3  4  6    (+) b b b b a b
  25       1  2  3  5  6    4                (+) a a a b a a
  26       4                1  2  3  5  6    (-) b b b a b b
  27       1  2  4  5  6    3                (-) a a b a a a
  28       3                1  2  3  5  6    (+) b b a b b b
  29       1  3  4  5  6    2                (+) a b a a a a
  30       2                1  3  4  5  6    (-) b a b b b b
  31       2  3  4  5  6    1                (-) b a a a a a
  32       1                2  3  4  5  6    (+) a b b b b b
*/
  
  if ( ioutfg.ndeti != 32 ) {
    cerr << "Inconsistent ndeti != 32 " << ioutfg.ndeti << " Aborting " << endl;
    GA::Terminate();
  } 
  
  if ( ioutfg.singly.size() != 6 ) {
    cerr << "Inconsistent number of open shells" << ioutfg.singly.size();
    cerr << " configuration index " << ioutfg.index  << " Aborting " << endl;
    GA::Terminate();
  } 
  
  const int pad = 2;
  vector<int> jdbl;
  jdbl.reserve( nalpha + pad );
  
  vector< vector<int> > ifiga( ndeti, vector<int>(nalpha) );
  vector< vector<int> > ifigb( ndeti, vector<int>(nalpha) );
  vector<int> ms( ndeti );
  
  vector<int>::iterator dit; 
  
  for (dit = ioutfg.doubly.begin(); dit != ioutfg.doubly.end(); ++dit ) {
    jdbl.push_back( (*dit) );
  }
  
  // determinants with ms = 0
  
  int ndeta = 1;
  int ndetb, ndetc;
  int num_doubly = ioutfg.doubly.size();
  
  for(int i3=3; i3<=6; ++i3) {
    jdbl[nalpha-1] = ioutfg.singly[i3-1];
    for( int i2=2; i2<=i3-1; ++i2) {
      jdbl[nalpha-2] = ioutfg.singly[i2-1];
      for( int i1=1; i1<=i2-1; ++i1) {
	jdbl[nalpha-3] = ioutfg.singly[i1-1];
        
	ndetb = min( ndeta, 41-ndeta );
	ndetc = min( ndeta+1, 40-ndeta);
	ms[ndetb-1] = msbas;
	
	//ifiga[ndetb-1].resize(nalpha);
	//ifigb[ndetc-1].resize(nalpha);
	for( int jj=1; jj<= nalpha; ++jj ) {
	  ifiga[ndetb-1][jj-1] = jdbl[jj-1];
	  ifigb[ndetc-1][jj-1] = jdbl[jj-1];
	}
	ndeta += 2;
      }
    }
  }
  
  // determinants with ms = 2, -2
  
  nalpha += 2;
  ndeta = 21;
  
  for(int i5=5; i5<=6; ++i5) {
    jdbl[nalpha-1] = ioutfg.singly[i5-1];
    for(int i4=4; i4<=i5-1; ++i4) {
      jdbl[nalpha-2] = ioutfg.singly[i4-1];
      for(int i3=3;i3<=i4-1; ++i3) {
        jdbl[nalpha-3] = ioutfg.singly[i3-1];
        for(int i2=2;i2<=i3-1; ++i2) {
          jdbl[nalpha-4] = ioutfg.singly[i2-1];
          for(int i1=1;i1<=i2-1; ++i1) {
            jdbl[nalpha-5] = ioutfg.singly[i1-1];
            ms[ndeta-1] = msbas + 2;
            ms[ndeta] = msbas - 2;
	    
            ifiga[ndeta-1].resize( nalpha );
            ifigb[ndeta].resize( nalpha );

            for (int jj=1; jj<=nalpha; ++jj) {
	      ifiga[ndeta-1][jj-1] = jdbl[jj-1];
	      ifigb[ndeta][jj-1] = jdbl[jj-1];
            }
            ndeta += 2;
          }
        }
      }
    }
  }
  
  nalpha -= 4;
  ndeta = 32;
  
  for(int i1=1; i1<=6; ++i1) {
    jdbl[nalpha-1] = ioutfg.singly[i1-1];
    
  ifiga[ndeta-1].resize(nalpha);
  ifigb[ndeta-2].resize(nalpha);

    for( int jj=1; jj<= nalpha; ++jj ) {
      ifiga[ndeta-1][jj-1] = jdbl[jj-1];
      ifigb[ndeta-2][jj-1] = jdbl[jj-1];
    }
    ndeta -= 2;
  }
  
  if ( pushBackToJoutfg( ioutfg, ifiga, ifigb, ms ) != 0 ) {
    cerr << "Error: fatal error from pushBackToJoutfg: << Aborting " << endl;
    GA::Terminate();
  }
  return( 0 );
}

//Eight open shell electrons. EVEN symmetry

int PsociDeterminants::eightOpenShellEvenDeterminants( JOUTFG & ioutfg )
{ 
  int msbas = fetchMaxOpenShells() / 2;
  int nalpha = l_nelec / 2;
  int ndeti = ioutfg.ndeti;
  
/*
  case 5:  eight open-shell electrons

          number of double-group functions =  64
          number of determinants           = 128

  arrangement of the singly occupied orbitals is analogous to
  that for case 4.
*/
    
  if ( ioutfg.ndeti != 128 ) {
    cerr << "Inconsistent ndeti != 128 " << ioutfg.ndeti << " Aborting " << endl;
    GA::Terminate();
  }
  
  if ( ioutfg.singly.size() != 8 ) {
    cerr << "Inconsistent number of open shells" << ioutfg.singly.size();
    cerr << " configuration index " << ioutfg.index  << " Aborting " << endl;
    GA::Terminate();
  }
  
  const int pad = 4;
  vector<int> jdbl;
  jdbl.reserve( nalpha + pad );
  
  vector< vector<int> > ifiga( ndeti, vector<int>(nalpha) );
  vector< vector<int> > ifigb( ndeti, vector<int>(nalpha) );
  vector<int> ms( ndeti );
  
  vector<int>::iterator dit;
  
  for (dit = ioutfg.doubly.begin(); dit != ioutfg.doubly.end(); ++dit ) {
    jdbl.push_back( (*dit) );
  }
  
  int ndeta = 1;
  int ndetb, ndetc;
  int num_doubly = ioutfg.doubly.size();
  
  // determinants with ms = 0
  
  for(int i4=4; i4<=8; ++i4) {
    jdbl[nalpha-1] = ioutfg.singly[i4-1];
    for(int i3=3; i3<=i4-1; ++i3) {
      jdbl[nalpha-2] = ioutfg.singly[i3-1];
      for( int i2=2; i2<=i3-1; ++i2) {
	jdbl[nalpha-3] = ioutfg.singly[i2-1];
	for( int i1=1; i1<=i2-1; ++i1) {
	  jdbl[nalpha-4] = ioutfg.singly[i1-1];
	  
	  ndetb = min( ndeta, 141-ndeta );
	  ndetc = min( ndeta+1, 140-ndeta);
	  ms[ndetb-1] = msbas;
	  
	  //ifiga[ndetb-1].resize(nalpha);
	  //ifigb[ndetc-1].resize(nalpha);
	  for( int jj=1; jj<= nalpha; ++jj ) {
	    ifiga[ndetb-1][jj-1] = jdbl[jj-1];
	    ifigb[ndetc-1][jj-1] = jdbl[jj-1];
	  }
	  ndeta += 2;
	}
      }
    }
  }
  
  
  // determinants with ms = 2, -2
  
  nalpha += 2;
  ndeta = 71;
  
  for(int i6=6; i6<=8; ++i6) {
    jdbl[nalpha-1] = ioutfg.singly[i6-1];
    for(int i5=5; i5<=i6-1; ++i5) {
      jdbl[nalpha-2] = ioutfg.singly[i5-1]; 
      for(int i4=4; i4<=i5-1; ++i4) {
	jdbl[nalpha-3] = ioutfg.singly[i4-1];
	for(int i3=3;i3<=i4-1; ++i3) {
	  jdbl[nalpha-4] = ioutfg.singly[i3-1];
	  for(int i2=2;i2<=i3-1; ++i2) {
	    jdbl[nalpha-5] = ioutfg.singly[i2-1];
	    for(int i1=1;i1<=i2-1; ++i1) {
	      jdbl[nalpha-6] = ioutfg.singly[i1-1];
	      ms[ndeta-1] = msbas + 2;
	      ms[ndeta] = msbas - 2;
	      
	      ifiga[ndeta-1].resize( nalpha );
	      ifigb[ndeta].resize( nalpha );
	      
	      for (int jj=1; jj<=nalpha; ++jj) {
		ifiga[ndeta-1][jj-1] = jdbl[jj-1];
		ifigb[ndeta][jj-1] = jdbl[jj-1];
	      }
	      ndeta += 2;
	    }
	  }
	}
      }
    }
  }
  
  nalpha -= 4;
  ndeta = 126;
  
  for(int i2=2; i2<=8; ++i2) {
    jdbl[nalpha-1] = ioutfg.singly[i2-1];
    for(int i1=1; i1<=i2-1; ++i1) {
      jdbl[nalpha-2] = ioutfg.singly[i1-1];
      
      ifiga[ndeta-1].resize(nalpha);
      ifigb[ndeta-2].resize(nalpha);
      
      for( int jj=1; jj<= nalpha; ++jj ) {
	ifiga[ndeta-1][jj-1] = jdbl[jj-1];
	ifigb[ndeta-2][jj-1] = jdbl[jj-1];
      }
      ndeta -= 2; 
    }
  }
  
  // determinants with ms = 4, -4 
  
  ms[126] = msbas + 4;
  ms[127] = msbas - 4;
  
  for(int i1=1; i1<=8; ++i1) {
    jdbl[num_doubly+i1-1] = ioutfg.singly[i1-1];
  }

  ifiga[127].resize(num_doubly);
  ifigb[126].resize(num_doubly);
  for(int jj=1; jj<=num_doubly; ++jj) {
    ifiga[127][jj-1] = jdbl[jj-1];
    ifigb[126][jj-1] = jdbl[jj-1];
  }

  nalpha += 6;
  
  ifiga[126].resize( nalpha );
  ifigb[127].resize( nalpha );
  for(int jj=1; jj<=nalpha; ++jj) {
    ifigb[127][jj-1] = jdbl[jj-1];
    ifiga[126][jj-1] = jdbl[jj-1];
  }
  
  if ( pushBackToJoutfg( ioutfg, ifiga, ifigb, ms ) != 0 ) {
    cerr << "Error: fatal error from pushBackToJoutfg: << Aborting " << endl;
    GA::Terminate();
  }
  return( 0 );
}

//10 open shell electrons EVEN symmetry

int PsociDeterminants::tenOpenShellEvenDeterminants( JOUTFG & ioutfg )
{ 
  int msbas = fetchMaxOpenShells() / 2;
  int nalpha = l_nelec / 2;
  int ndeti = ioutfg.ndeti;
  
/*
  case 6:  ten open-shell electrons

          number of double-group functions = 256
          number of determinants           = 512

  arrangement of the singly occupied orbitals is analogous to
  that for case 4.
*/
    
  if ( ioutfg.ndeti != 512 ) {
    cerr << "Inconsistent ndeti != 512 " << ioutfg.ndeti << " Aborting " << endl;
    GA::Terminate();
  }
  
  if ( ioutfg.singly.size() != 10 ) {
    cerr << "Inconsistent number of open shells" << ioutfg.singly.size();
    cerr << " configuration index " << ioutfg.index  << " Aborting " << endl;
    GA::Terminate();
  }
  
  const int pad = 4;
  vector<int> jdbl;
  jdbl.reserve( nalpha + pad );
  
  vector< vector<int> > ifiga( ndeti, vector<int>(nalpha) );
  vector< vector<int> > ifigb( ndeti, vector<int>(nalpha) );
  vector<int> ms( ndeti );
  
  vector<int>::iterator dit;
  
  for (dit = ioutfg.doubly.begin(); dit != ioutfg.doubly.end(); ++dit ) {
    jdbl.push_back( (*dit) );
  }
  
  int ndeta = 1;
  int ndetb, ndetc;
  int num_doubly = ioutfg.doubly.size();
  
  // determinants with ms = 0
  
  for(int i5=5; i5<=10; ++i5) {
    jdbl[nalpha-1] = ioutfg.singly[i5-1];
    for(int i4=4; i4<=i5-1; ++i4) {
      jdbl[nalpha-2] = ioutfg.singly[i4-1];
      for(int i3=3; i3<=i4-1; ++i3) {
	jdbl[nalpha-3] = ioutfg.singly[i3-1];
	for( int i2=2; i2<=i3-1; ++i2) {
	  jdbl[nalpha-4] = ioutfg.singly[i2-1];
	  for( int i1=1; i1<=i2-1; ++i1) {
	    jdbl[nalpha-5] = ioutfg.singly[i1-1];
	    
	    ndetb = min( ndeta, 505-ndeta );
	    ndetc = min( ndeta+1, 504-ndeta);
	    ms[ndetb-1] = msbas;
	    
	    //ifiga[ndetb-1].resize(nalpha);
	    //ifigb[ndetc-1].resize(nalpha);
	    for( int jj=1; jj<= nalpha; ++jj ) {
	      ifiga[ndetb-1][jj-1] = jdbl[jj-1];
	      ifigb[ndetc-1][jj-1] = jdbl[jj-1];
	    }
	    ndeta += 2;
	  }
	}
      }
    }
  }
  
 // determinants with ms = 2, -2
  
  nalpha += 2;
  ndeta = 253;
  
  for(int i7=7; i7<=10; ++i7) {
    jdbl[nalpha-1] = ioutfg.singly[i7-1];
    for(int i6=6; i6<=i7-1; ++i6) {
      jdbl[nalpha-2] = ioutfg.singly[i6-1];
      for(int i5=5; i5<=i6-1; ++i5) {
	jdbl[nalpha-3] = ioutfg.singly[i5-1]; 
	for(int i4=4; i4<=i5-1; ++i4) {
	  jdbl[nalpha-4] = ioutfg.singly[i4-1];
	  for(int i3=3;i3<=i4-1; ++i3) {
	    jdbl[nalpha-5] = ioutfg.singly[i3-1];
	    for(int i2=2;i2<=i3-1; ++i2) {
	      jdbl[nalpha-6] = ioutfg.singly[i2-1];
	      for(int i1=1;i1<=i2-1; ++i1) {
		jdbl[nalpha-7] = ioutfg.singly[i1-1];
		
		ms[ndeta-1] = msbas + 2;
		ms[ndeta] = msbas - 2;
		
		ifiga[ndeta-1].resize( nalpha );
		ifigb[ndeta].resize( nalpha );
		
		for (int jj=1; jj<=nalpha; ++jj) {
		  ifiga[ndeta-1][jj-1] = jdbl[jj-1];
		  ifigb[ndeta][jj-1] = jdbl[jj-1];
		}
		ndeta += 2;
	      }
	    }
	  }
	}
      }
    }
  }
  
  nalpha -= 4;
  ndeta = 492;
  
  for(int i3=3; i3<=10; ++i3) {
    jdbl[nalpha-1] = ioutfg.singly[i3-1];
    for(int i2=2; i2<=i3-1; ++i2) {
      jdbl[nalpha-2] = ioutfg.singly[i2-1];
      for(int i1=1; i1<=i2-1; ++i1) {
	jdbl[nalpha-3] = ioutfg.singly[i1-1];
	
	ifiga[ndeta-1].resize(nalpha);
	ifigb[ndeta-2].resize(nalpha);
	
	for( int jj=1; jj<= nalpha; ++jj ) {
	  ifiga[ndeta-1][jj-1] = jdbl[jj-1];
	  ifigb[ndeta-2][jj-1] = jdbl[jj-1];
	}
	ndeta -= 2;
      }
    }
  }
  
// determinants with ms = 4, -4

  nalpha += 6;
  ndeta = 493;
  
  for(int i9=9; i9<=10; ++i9) {
    jdbl[nalpha-1] = ioutfg.singly[i9-1];
    for(int i8=8; i8<=i9-1; ++i8) {
      jdbl[nalpha-2] = ioutfg.singly[i8-1];
      for(int i7=7; i7<=i8-1; ++i7) {
	jdbl[nalpha-3] = ioutfg.singly[i7-1];
	for(int i6=6; i6<=i7-1; ++i6) {
	  jdbl[nalpha-4] = ioutfg.singly[i6-1];
	  for(int i5=5; i5<=i6-1; ++i5) {
	    jdbl[nalpha-5] = ioutfg.singly[i5-1]; 
	    for(int i4=4; i4<=i5-1; ++i4) {
	      jdbl[nalpha-6] = ioutfg.singly[i4-1];
	      for(int i3=3;i3<=i4-1; ++i3) {
		jdbl[nalpha-7] = ioutfg.singly[i3-1];
		for(int i2=2;i2<=i3-1; ++i2) {
		  jdbl[nalpha-8] = ioutfg.singly[i2-1];
		  for(int i1=1;i1<=i2-1; ++i1) {
		    jdbl[nalpha-9] = ioutfg.singly[i1-1];
		    
		    ms[ndeta-1] = msbas + 4;
		    ms[ndeta] = msbas - 4;
		    
		    ifiga[ndeta-1].resize( nalpha );
		    ifigb[ndeta].resize( nalpha );
		    
		    for (int jj=1; jj<=nalpha; ++jj) {
		      ifiga[ndeta-1][jj-1] = jdbl[jj-1];
		      ifigb[ndeta][jj-1] = jdbl[jj-1];
		    }
		    ndeta += 2;
		  }
		}
	      }
	    }
	  }
	}
      }
    }
  }
  
  nalpha -= 8;
  ndeta = 512;
  
  for(int i1=1; i1<=10; ++i1) {
    jdbl[nalpha-1] = ioutfg.singly[i1-1];
    
    ifiga[ndeta-1].resize(nalpha);
    ifigb[ndeta-2].resize(nalpha);
    
    for( int jj=1; jj<= nalpha; ++jj ) {
      ifiga[ndeta-1][jj-1] = jdbl[jj-1];
      ifigb[ndeta-2][jj-1] = jdbl[jj-1];
    }
    ndeta -= 2; 
  }
  
  if ( pushBackToJoutfg( ioutfg, ifiga, ifigb, ms ) != 0 ) {
    cerr << "Error: fatal error from pushBackToJoutfg: << Aborting " << endl;
    GA::Terminate();
  }
  return( 0 );
}

// ODD symmetry determinants
/*
  generate all of the determinants with ms odd needed to make
  double-group-adapted functions from a spatial configuratio
*/

/* TODO: The original code did a kludgy thing for this piece. It passed the 
   building to the twoOpenShellOdd code and passed empty orbitals to dets 1 and 2.
   I have no way to confirm why this was done. It may need to be revisited. For now
   simply construct a single det, with alpha/beta == doubly_occ orbs

*/

// Formally called MSODET
int PsociDeterminants::zeroOpenShellOddDeterminants( JOUTFG & ioutfg )
{
  int msbas = fetchMaxOpenShells() / 2; 
  // int nalpha = 1 + ( l_nelec ) / 2; 
  int ndeti = ioutfg.ndeti;
  
  /*
    case 1:  zero open-shell electron
    
    number of double-group functions = 1
    number of determinants           = 1
    
    the list of doubly occupied orbitals goes in both the spin-up
    and spin-down lists. the singly occupied orbital goes in the
    spin-up list.
    
  */
  
  if ( ioutfg.ndeti != 1 ) {
    cerr << "Inconsistent ndeti != 1 " << ioutfg.ndeti << " Aborting " << endl;
    GA::Terminate();
  }
  
  if ( ioutfg.singly.size() != 0 ) {
    cerr << "Inconsistent number of open shells" << ioutfg.singly.size();
    cerr << " configuration index " << ioutfg.index  << " Aborting " << endl;
    GA::Terminate();
  }
  
  int num_doubly = ioutfg.doubly.size();
  
  vector<int> jdbl;
  jdbl.reserve( num_doubly );
  
  vector< vector<int> > ifiga( ndeti, vector<int>(num_doubly) );
  vector< vector<int> > ifigb( ndeti, vector<int>(num_doubly) );
  vector<int> ms( ndeti );
  
  // Copy all current double occupied orbs
  vector<int>::iterator dit;
  
  for (dit = ioutfg.doubly.begin(); dit != ioutfg.doubly.end(); ++dit ) {
    jdbl.push_back( (*dit) );
  }
  
  int ndeta = 1;
  ms[0] = msbas;
  
  ifiga[ndeta-1].resize(num_doubly);
  ifigb[ndeta-1].resize(num_doubly);
  for( int jj=1; jj<= num_doubly; ++jj ) {
    ifiga[ndeta-1][jj-1] = jdbl[jj-1];
    ifigb[ndeta-1][jj-1] = jdbl[jj-1];
  }
  
  if ( pushBackToJoutfg( ioutfg, ifiga, ifigb, ms ) != 0 ) {
    cerr << "Error: fatal error from pushBackToJoutfg: << Aborting " << endl;
    GA::Terminate();
  }
  return( 0 );
}   
int PsociDeterminants::twoOpenShellOddDeterminants( JOUTFG & ioutfg )
{
  int msbas = fetchMaxOpenShells() / 2; 
  int nalpha = 1 + l_nelec / 2; 
  int ndeti = ioutfg.ndeti;
  
  /*
    case 1:  two open-shell electrons
    
    number of double-group functions = 1
    number of determinants           = 2
    
    arrangement of the singly occupied orbitals
    
    singly occupied orbitals
    det   alpha list   beta list    (sign) spin function
    
    1       1  2                    (+) a a
    2                    1  2       (+) b b
  */
  
  if ( ioutfg.ndeti != 2 ) {
    cerr << "Inconsistent ndeti != 2 " << ioutfg.ndeti << " Aborting " << endl;
    GA::Terminate();
  }
  
  if ( ioutfg.singly.size() != 2 ) {
    cerr << "Inconsistent number of open shells" << ioutfg.singly.size();
    cerr << " configuration index " << ioutfg.index  << " Aborting " << endl;
    GA::Terminate();
  }
  
  vector<int> jdbl;
  jdbl.reserve( nalpha );
  
  vector< vector<int> > ifiga( ndeti, vector<int>(nalpha) );
  vector< vector<int> > ifigb( ndeti, vector<int>(nalpha) );
  vector<int> ms( ndeti );
  
  vector<int>::iterator dit;
  
  for (dit = ioutfg.doubly.begin(); dit != ioutfg.doubly.end(); ++dit ) {
    jdbl.push_back( (*dit) );
  }
  
  int ndeta = 1;
  int num_doubly = ioutfg.doubly.size();

  ms[0] = msbas + 1;
  ms[1] = msbas - 1;

  jdbl[num_doubly] = ioutfg.singly[0];
  jdbl[num_doubly + 1] = ioutfg.singly[1];
  
  ifiga[ndeta].resize(num_doubly); 
  ifigb[ndeta-1].resize(num_doubly);
  for( int jj=1; jj<= num_doubly; ++jj ) {
    ifiga[1][jj-1] = jdbl[jj-1];
    ifigb[0][jj-1] = jdbl[jj-1];
  }
  for( int jj=1; jj<= nalpha; ++jj ) {
    ifiga[0][jj-1] = jdbl[jj-1];
    ifigb[1][jj-1] = jdbl[jj-1];
  }
  
  if ( pushBackToJoutfg( ioutfg, ifiga, ifigb, ms ) != 0 ) {
    cerr << "Error: fatal error from pushBackToJoutfg: << Aborting " << endl;
    GA::Terminate();
  }
  return( 0 );
}   

// Odd determinants
int PsociDeterminants::fourOpenShellOddDeterminants( JOUTFG & ioutfg )
{
  int msbas = fetchMaxOpenShells() / 2; 
  int nalpha = 1 + l_nelec / 2; 
  int ndeti = ioutfg.ndeti;
  
  /*
    number of double-group functions = 4
    number of determinants           = 8
    
    arrangement of the singly occupied orbitals
    
    singly occupied orbitals
    det   alpha list   beta list    (sign) spin function
    
    1      1  2  3     4            (+) a a a b
    2      4           1  2  3      (-) b b b a
    3      1  2  4     3            (-) a a b a
    4      3           1  2  4      (+) b b a b
    5      1  3  4     2            (+) a b a a
    6      2           1  3  4      (-) b a b b
    7      2  3  4     1            (-) b a a a
    8      1           2  3  4      (+) a b b b
  */
  
  if ( ioutfg.ndeti != 8 ) {
    cerr << "Inconsistent ndeti != 8 " << ioutfg.ndeti << " Aborting " << endl;
    GA::Terminate();
  }
  
  if ( ioutfg.singly.size() != 4 ) {
    cerr << "Inconsistent number of open shells" << ioutfg.singly.size();
    cerr << " configuration index " << ioutfg.index  << " Aborting " << endl;
    GA::Terminate();
  }
  
  vector<int> jdbl;
  jdbl.reserve( nalpha );
  
  vector< vector<int> > ifiga( ndeti, vector<int>(nalpha) );
  vector< vector<int> > ifigb( ndeti, vector<int>(nalpha) );
  vector<int> ms( ndeti );
  
  int ndeta = 1;
  vector<int>::iterator dit;
  
  for (dit = ioutfg.doubly.begin(); dit != ioutfg.doubly.end(); ++dit ) {
    jdbl.push_back( (*dit) );
  }
  
  for( int i3=3; i3<=4; ++i3) {
    jdbl[nalpha-1] = ioutfg.singly[i3-1];
    for( int i2=2; i2<=i3-1; ++i2) {
      jdbl[nalpha-2] = ioutfg.singly[i2-1];
      for( int i1=1; i1<=i2-1; ++i1) {
	jdbl[nalpha-3] = ioutfg.singly[i1-1];
	
	ms[ndeta-1] = msbas + 1;
	ms[ndeta] = msbas - 1;
	
	for( int jj=1; jj<= nalpha; ++jj ) {
	  ifiga[ndeta-1][jj-1] = jdbl[jj-1];
	  ifigb[ndeta][jj-1] = jdbl[jj-1];
	}
	ndeta += 2;
      }
    }
  }
  
  nalpha -= 2;
  ndeta = 8;  
  
  for (int i1=1; i1<=4; ++i1) {
  ifiga[ndeta-1].resize(nalpha); 
  ifigb[ndeta-2].resize(nalpha);
    jdbl[nalpha-1] = ioutfg.singly[i1-1];
    for( int jj=1; jj<= nalpha; ++jj ) {
      ifiga[ndeta-1][jj-1] = jdbl[jj-1];
      ifigb[ndeta-2][jj-1] = jdbl[jj-1];
    }
    ndeta -= 2;
  }
  
  if ( pushBackToJoutfg( ioutfg, ifiga, ifigb, ms ) != 0 ) {
    cerr << "Error: fatal error from pushBackToJoutfg: << Aborting " << endl;
    GA::Terminate();
  }
  return( 0 );
}   

//Odd determinants
int PsociDeterminants::sixOpenShellOddDeterminants( JOUTFG & ioutfg )
{
  int msbas = fetchMaxOpenShells() / 2; 
  int nalpha = 1 + l_nelec / 2; 
  int ndeti = ioutfg.ndeti;
  
  /*
    number of double-group functions = 16
    number of determinants           = 32
    
    arrangement of the singly occupied orbitals
    
    singly occupied orbitals
    det       alpha list        beta list     (sign) spin function
    
    1     1  2  3  4       5  6               (+) a a a a b b
    2     5  6             1  2  3  4         (+) b b a a a a
    3     1  2  3  5       4  6               (-) a a a b a b
    4     4  6             1  2  3  5         (-) b b b a b a
    5     1  2  4  5       3  6               (+) a a b a a b
    6     3  6             1  2  4  5         (+) b b a b b a
    7     1  3  4  5       2  6               (-) a b a a a b
    8     2  6             1  3  4  5         (-) b a b b b a
    9     2  3  4  5       1  6               (+) b a a a a b
    10     1  6             2  3  4  5         (+) a b b b b a
    11     1  2  3  6       4  5               (+) a a a b b a
    12     4  5             1  2  3  6         (+) b b b a a b
    13     1  2  4  6       3  5               (-) a a b a b a
    14     3  5             1  2  4  6         (-) b b a b a b
    15     1  3  4  6       2  5               (+) a b a a b a
    16     2  5             1  3  4  6         (+) b a b b a b
    17     2  3  4  6       1  5               (-) b a a a b a
    18     1  5             2  3  4  6         (-) a b b b a b
    19     1  2  5  6       3  4               (+) a a b b a a
    20     3  4             1  2  5  6         (+) b b a a b b
    21     1  3  5  6       2  4               (-) a b a b a a
    22     2  4             1  3  5  6         (-) b a b a b b
    23     2  3  5  6       1  4               (+) b a a b a a
    24     1  4             2  3  5  6         (+) a b b a b b
    25     1  4  5  6       2  3               (+) a b b a a a
    26     2  3             1  4  5  6         (+) b a a b b b
    27     2  4  5  6       1  3               (-) b a b a a a
    28     1  3             2  4  5  6         (-) a b a b b b
    29     3  4  5  6       1  2               (+) b b a a a a
    30     1  2             3  4  5  6         (+) a a b b b b
    31     1  2  3  4  5  6                    (+) a a a a a a
    32                      1  2  3  4  5  6   (+) b b b b b b
  */
  
  if ( ioutfg.ndeti != 32 ) {
    cerr << "Inconsistent ndeti != 32 " << ioutfg.ndeti << " Aborting " << endl;
    GA::Terminate();
  }
  
  if ( ioutfg.singly.size() != 6 ) {
    cerr << "Inconsistent number of open shells" << ioutfg.singly.size();
    cerr << " configuration index " << ioutfg.index  << " Aborting " << endl;
    GA::Terminate();
  }
  
  const int pad = 2;
  vector<int> jdbl;
  jdbl.reserve( nalpha + pad );
  
  vector< vector<int> > ifiga( ndeti, vector<int>(nalpha) );
  vector< vector<int> > ifigb( ndeti, vector<int>(nalpha) );
  vector<int> ms( ndeti );
  
  int ndeta = 1;
  int num_doubly = ioutfg.doubly.size();

  vector<int>::iterator dit;
  
  for (dit = ioutfg.doubly.begin(); dit != ioutfg.doubly.end(); ++dit ) {
    jdbl.push_back( (*dit) );
  }
  
  //  determinants with ms = 1,-1
  
  for( int i4=4; i4<=6; ++i4) {
    jdbl[nalpha-1] = ioutfg.singly[i4-1];
    for( int i3=3; i3<=i4-1; ++i3) {
      jdbl[nalpha-2] = ioutfg.singly[i3-1];
      for( int i2=2; i2<=i3-1; ++i2) {
	jdbl[nalpha-3] = ioutfg.singly[i2-1];
	for( int i1=1; i1<=i2-1; ++i1) {
	  jdbl[nalpha-4] = ioutfg.singly[i1-1];
	  
	  ms[ndeta-1] = msbas + 1;
	  ms[ndeta] = msbas - 1;
	  
	  for( int jj=1; jj<= nalpha; ++jj ) {
	    ifiga[ndeta-1][jj-1] = jdbl[jj-1];
	    ifigb[ndeta][jj-1] = jdbl[jj-1];
	  }
	  ndeta += 2;
	}
      }
    }
  }
  
  nalpha -= 2;
  ndeta = 30;
  
  for( int i2=2; i2<=6; ++i2) {
    jdbl[nalpha-1] = ioutfg.singly[i2-1];
    for( int i1=1; i1<=i2-1; ++i1) {
      jdbl[nalpha-2] = ioutfg.singly[i1-1];
      
      ifiga[ndeta-1].resize(nalpha);
      ifigb[ndeta-2].resize(nalpha);
      
      for( int jj=1; jj<= nalpha; ++jj ) {
	ifiga[ndeta-1][jj-1] = jdbl[jj-1];
	ifigb[ndeta-2][jj-1] = jdbl[jj-1];
      }
      ndeta -= 2;
    }
  }
  
  // determinants with ms = 3,-3
  
  ms[30] = msbas + 3;
  ms[31] = msbas - 3;
  
  ifiga[31].resize(num_doubly);
  ifigb[30].resize(num_doubly);
  for(int jj=1; jj<=num_doubly; ++jj) {
    ifiga[31][jj-1] = jdbl[jj-1];
    ifigb[30][jj-1] = jdbl[jj-1];
  }
  for(int i1=1; i1<=6; ++i1) {
    jdbl[num_doubly+i1-1] = ioutfg.singly[i1-1];
  }
  
  nalpha += 4;
  
  ifiga[30].resize( nalpha );
  ifigb[31].resize( nalpha );
  for(int jj=1; jj<=nalpha; ++jj) {
    ifigb[31][jj-1] = jdbl[jj-1];
    ifiga[30][jj-1] = jdbl[jj-1];
  }
  
  if ( pushBackToJoutfg( ioutfg, ifiga, ifigb, ms ) != 0 ) {
    cerr << "Error: fatal error from pushBackToJoutfg: << Aborting " << endl;
    GA::Terminate();
  }
  return( 0 );
}   

//Odd determinants
int PsociDeterminants::eightOpenShellOddDeterminants( JOUTFG & ioutfg )
{
  int msbas = fetchMaxOpenShells() / 2; 
  int nalpha = 1 + l_nelec / 2; 
  int ndeti = ioutfg.ndeti;
  
  /*    
  number of double-group functions =  64
  number of determinants           = 128
  
  arrangement of the singly occupied orbitals is analogous to
  that for case 3.
  */
  
  if ( ioutfg.ndeti != 128 ) {
    cerr << "Inconsistent ndeti != 128 " << ioutfg.ndeti << " Aborting " << endl;
    GA::Terminate();
  }
  
  if ( ioutfg.singly.size() != 8 ) {
    cerr << "Inconsistent number of open shells" << ioutfg.singly.size();
    cerr << " configuration index " << ioutfg.index  << " Aborting " << endl;
    GA::Terminate();
  }
  
  const int pad = 2;
  vector<int> jdbl;
  jdbl.reserve( nalpha + pad );
  
  vector< vector<int> > ifiga( ndeti, vector<int>(nalpha) );
  vector< vector<int> > ifigb( ndeti, vector<int>(nalpha) );
  vector<int> ms( ndeti );
  
  int ndeta = 1;
  vector<int>::iterator dit;
  
  for (dit = ioutfg.doubly.begin(); dit != ioutfg.doubly.end(); ++dit ) {
    jdbl.push_back( (*dit) );
  }
  
  //  determinants with ms = 1,-1
  
  for( int i5=5; i5<=8; ++i5) {
    jdbl[nalpha-1] = ioutfg.singly[i5-1];
    for( int i4=4; i4<=i5-1; ++i4) {
      jdbl[nalpha-2] = ioutfg.singly[i4-1];
      for( int i3=3; i3<=i4-1; ++i3) {
	jdbl[nalpha-3] = ioutfg.singly[i3-1];
	for( int i2=2; i2<=i3-1; ++i2) {
	  jdbl[nalpha-4] = ioutfg.singly[i2-1];
	  for( int i1=1; i1<=i2-1; ++i1) {
	    jdbl[nalpha-5] = ioutfg.singly[i1-1];
	    
	    ms[ndeta-1] = msbas + 1;
	    ms[ndeta] = msbas - 1;
	    
	    for( int jj=1; jj<= nalpha; ++jj ) {
	      ifiga[ndeta-1][jj-1] = jdbl[jj-1];
	      ifigb[ndeta][jj-1] = jdbl[jj-1];
	    }
	    ndeta += 2;
	  }
	}
      }
    }
  }
  
  nalpha -= 2;
  ndeta = 112;
  
  for( int i3=3; i3<=8; ++i3) {
    jdbl[nalpha-1] = ioutfg.singly[i3-1];
    for( int i2=2; i2<=i3-1; ++i2) {
      jdbl[nalpha-2] = ioutfg.singly[i2-1];
      for( int i1=1; i1<=i2-1; ++i1) {
	jdbl[nalpha-3] = ioutfg.singly[i1-1];
	
	ifiga[ndeta-1].resize( nalpha );
	ifigb[ndeta-2].resize( nalpha );
	
	for( int jj=1; jj<= nalpha; ++jj ) {
	  ifiga[ndeta-1][jj-1] = jdbl[jj-1];
	  ifigb[ndeta-2][jj-1] = jdbl[jj-1];
	}
	ndeta -= 2;
      }
    }
  }
  
  // determinants with ms = 3,-3
  
  nalpha += 4;
  ndeta = 113;
  
  for(int i7=7; i7<=8; ++i7) {
    jdbl[nalpha-1] = ioutfg.singly[i7-1];
    for(int i6=6; i6<=i7-1; ++i6) {
      jdbl[nalpha-2] = ioutfg.singly[i6-1];
      for(int i5=5; i5<=i6-1; ++i5) {
	jdbl[nalpha-3] = ioutfg.singly[i5-1];
	for(int i4=4; i4<=i5-1; ++i4) {
	  jdbl[nalpha-4] = ioutfg.singly[i4-1];
	  for(int i3=3;i3<=i4-1; ++i3) {
	    jdbl[nalpha-5] = ioutfg.singly[i3-1];
	    for(int i2=2;i2<=i3-1; ++i2) {
	      jdbl[nalpha-6] = ioutfg.singly[i2-1];
	      for(int i1=1;i1<=i2-1; ++i1) {
		jdbl[nalpha-7] = ioutfg.singly[i1-1];
		ms[ndeta-1] = msbas + 3;
		ms[ndeta] = msbas - 3;
		
		ifiga[ndeta-1].resize( nalpha );
		ifigb[ndeta].resize( nalpha );
		
		for (int jj=1; jj<=nalpha; ++jj) {
		  ifiga[ndeta-1][jj-1] = jdbl[jj-1];
		  ifigb[ndeta][jj-1] = jdbl[jj-1];
		}
		ndeta += 2;
	      }
	    }
	  }
	}
      }
    }
  }
  
  nalpha -= 6;
  ndeta = 128;
  
  for(int i1=1; i1<=8; ++i1) {
    ifiga[ndeta-1].resize(nalpha);
    ifigb[ndeta-2].resize(nalpha);
    jdbl[nalpha-1] = ioutfg.singly[i1-1];
    for( int jj=1; jj<= nalpha; ++jj ) {
      ifiga[ndeta-1][jj-1] = jdbl[jj-1];
      ifigb[ndeta-2][jj-1] = jdbl[jj-1];
    }
    ndeta -= 2;
  }
  
  if ( pushBackToJoutfg( ioutfg, ifiga, ifigb, ms ) != 0 ) {
    cerr << "Error: fatal error from pushBackToJoutfg: << Aborting " << endl;
    GA::Terminate();
  }
  return( 0 );
}   

//Odd dets
int PsociDeterminants::tenOpenShellOddDeterminants( JOUTFG & ioutfg )
{
  int msbas = fetchMaxOpenShells() / 2; 
  int nalpha = 1 + l_nelec / 2; 
  int ndeti = ioutfg.ndeti;
  
  /*
    number of double-group functions = 256 
    number of determinants           = 512 
    
    arrangement of the singly occupied orbitals is analogous to
    that for case 3.
  */
  
  if ( ioutfg.ndeti != 512 ) {
    cerr << "Inconsistent ndeti != 512 " << ioutfg.ndeti << " Aborting " << endl;
    GA::Terminate();
  }
  
  if ( ioutfg.singly.size() != 10 ) {
    cerr << "Inconsistent number of open shells" << ioutfg.singly.size();
    cerr << " configuration index " << ioutfg.index  << " Aborting " << endl;
    GA::Terminate();
  }
  
  const int pad = 4;
  vector<int> jdbl;
  jdbl.reserve( nalpha + pad );
  
  vector< vector<int> > ifiga( ndeti, vector<int>(nalpha) );
  vector< vector<int> > ifigb( ndeti, vector<int>(nalpha) );
  vector<int> ms( ndeti );
  
  int ndeta = 1;
  int num_doubly = ioutfg.doubly.size();

  vector<int>::iterator dit;
  
  for (dit = ioutfg.doubly.begin(); dit != ioutfg.doubly.end(); ++dit ) {
    jdbl.push_back( (*dit) );
  }
  
  //  determinants with ms = 1,-1


  for( int i6=6; i6<=10; ++i6) {
    jdbl[nalpha-1] = ioutfg.singly[i6-1];
    for( int i5=5; i5<=i6-1; ++i5) {
      jdbl[nalpha-2] = ioutfg.singly[i5-1];
      for( int i4=4; i4<=i5-1; ++i4) {
	jdbl[nalpha-3] = ioutfg.singly[i4-1];
	for( int i3=3; i3<=i4-1; ++i3) {
	  jdbl[nalpha-4] = ioutfg.singly[i3-1];
	  for( int i2=2; i2<=i3-1; ++i2) {
	    jdbl[nalpha-5] = ioutfg.singly[i2-1];
	    for( int i1=1; i1<=i2-1; ++i1) {
	      jdbl[nalpha-6] = ioutfg.singly[i1-1];
	      
	      ms[ndeta-1] = msbas + 1;
	      ms[ndeta] = msbas - 1;
	      
	      for( int jj=1; jj<= nalpha; ++jj ) {
		ifiga[ndeta-1][jj-1] = jdbl[jj-1];
		ifigb[ndeta][jj-1] = jdbl[jj-1];
	      }
	      ndeta += 2;
	    }
	  }
	}
      }
    }
  }
  
  nalpha -= 2;
  ndeta = 420;
  
  for( int i4=4; i4<=10; ++i4) {
    jdbl[nalpha-1] = ioutfg.singly[i4-1];
    for( int i3=3; i3<=i4-1; ++i3) {
      jdbl[nalpha-2] = ioutfg.singly[i3-1];
      for( int i2=2; i2<=i3-1; ++i2) {
	jdbl[nalpha-3] = ioutfg.singly[i2-1];
	for( int i1=1; i1<=i2-1; ++i1) {
	  jdbl[nalpha-4] = ioutfg.singly[i1-1];
	  
	  ifiga[ndeta-1].resize( nalpha );
	  ifigb[ndeta-2].resize( nalpha );
	  
	  for( int jj=1; jj<= nalpha; ++jj ) {
	    ifiga[ndeta-1][jj-1] = jdbl[jj-1];
	    ifigb[ndeta-2][jj-1] = jdbl[jj-1];
	  }
	  ndeta -= 2;
	}
      }
    }
  }
  
  // determinants with ms = 3,-3
  
  nalpha += 4;
  ndeta = 421;
  
  for(int i8=8; i8<=10; ++i8) {
    jdbl[nalpha-1] = ioutfg.singly[i8-1];
    for(int i7=7; i7<=i8-1; ++i7) {
      jdbl[nalpha-2] = ioutfg.singly[i7-1];
      for(int i6=6; i6<=i7-1; ++i6) {
	jdbl[nalpha-3] = ioutfg.singly[i6-1];
	for(int i5=5; i5<=i6-1; ++i5) {
	  jdbl[nalpha-4] = ioutfg.singly[i5-1];
	  for(int i4=4; i4<=i5-1; ++i4) {
	    jdbl[nalpha-5] = ioutfg.singly[i4-1];
	    for(int i3=3;i3<=i4-1; ++i3) {
	      jdbl[nalpha-6] = ioutfg.singly[i3-1];
	      for(int i2=2;i2<=i3-1; ++i2) {
		jdbl[nalpha-7] = ioutfg.singly[i2-1];
		for(int i1=1;i1<=i2-1; ++i1) {
		  jdbl[nalpha-8] = ioutfg.singly[i1-1];
		  ms[ndeta-1] = msbas + 3;
		  ms[ndeta] = msbas - 3;
		  
		  ifiga[ndeta-1].resize( nalpha );
		  ifigb[ndeta].resize( nalpha );
		  
		  for (int jj=1; jj<=nalpha; ++jj) {
		    ifiga[ndeta-1][jj-1] = jdbl[jj-1];
		    ifigb[ndeta][jj-1] = jdbl[jj-1];
		  }
		  ndeta += 2;
		}
	      }
	    }
	  }
	}
      }
    }
  }
  
  nalpha -= 6;
  ndeta = 510;

  for(int i2=2; i2<=10; ++i2) {
    jdbl[nalpha-1] = ioutfg.singly[i2-1];
    for(int i1=1; i1<=i2-1; ++i1) {
      jdbl[nalpha-2] = ioutfg.singly[i1-1];
      ifiga[ndeta-1].resize(nalpha);
      ifigb[ndeta-2].resize(nalpha);
      
      for( int jj=1; jj<= nalpha; ++jj ) {
	ifiga[ndeta-1][jj-1] = jdbl[jj-1];
	ifigb[ndeta-2][jj-1] = jdbl[jj-1];
      }
      ndeta -= 2;
    }
  }
  
// determinants with ms = 5,-5

  ms[510] = msbas + 5;
  ms[511] = msbas - 5;
  
  ifiga[511].resize(num_doubly);
  ifigb[510].resize(num_doubly);
  for(int jj=1; jj<=num_doubly; ++jj) {
    ifiga[511][jj-1] = jdbl[jj-1];
    ifigb[510][jj-1] = jdbl[jj-1];
  }
  for(int i1=1; i1<=10; ++i1) {
    jdbl[num_doubly+i1-1] = ioutfg.singly[i1-1];
  }
  nalpha += 8;
  
  ifiga[510].resize( nalpha );
  ifigb[511].resize( nalpha );
  for(int jj=1; jj<=nalpha; ++jj) {
    ifigb[511][jj-1] = jdbl[jj-1];
    ifiga[510][jj-1] = jdbl[jj-1];
  }
  
  if ( pushBackToJoutfg( ioutfg, ifiga, ifigb, ms ) != 0 ) {
    cerr << "Error: fatal error from pushBackToJoutfg: << Aborting " << endl;
    GA::Terminate();
  }
  return( 0 );
}   

/* Compute the phase of the determinants by counting number of permutations required
   to bring to normal order. Then sort orbitals in ascending order

   Also compute sums: NOTE: sums must be relative to a Fortran indexing (1,2,,nbf)

   Treat all alpha/beta for a given iconfig
*/

int PsociDeterminants::computePhaseAndSort( JOUTFG & ioutfg ) {
  int nbf = l_nbf;
  int num_processed=0;
  
#ifdef DETAILEDCHECK
  if ( ioutfg.phases.size() != 0 ) ioutfg.phases.clear();
  if ( ioutfg.sumsa.size() != 0 ) ioutfg.sumsa.clear();
  if ( ioutfg.sumsb.size() != 0 ) ioutfg.sumsb.clear();
#endif
  
/* HPCToolkit suggests removing this
  if ( nbf <= 0 ) {
    cerr << "Error: unreasonable value for NBF was " << nbf << " Aborting " << endl;
    GA::Terminate();
  }
*/
  
  vector<vector<int> >::iterator ideter;
  vector<int>::iterator it;
  vector<int> phasea;
  vector<int> phaseb;
  
  // Alpha orbitals first
  
  /* NOTE: Orbital values are indexed (0,nbf-1). This impacts usage of lsum which presumes (1,nbf)
   */
  
/* HPCToolkit indicates too much time here
   Need to move to a more efficient algorithm
*/

  for( ideter=ioutfg.alpha_spins.begin(); ideter != ioutfg.alpha_spins.end(); ++ideter ) {
    
    int norbs = (*ideter).size();
    
    vector<int> buffer(nbf,0);
    int lsuma = 0; 
    int isign = 0;
    
    for(it = (*ideter).begin(); it != (*ideter).end(); ++it ) {
      int orb = (*it);
      lsuma += orb+1; 
      buffer[orb] = 1;
      ++num_processed;
      
      if ( orb != nbf-1 ) {
	for( int i=orb+1; i<nbf; ++i) {
	  isign += buffer[i];
	}
      }
    }
    ioutfg.sumsa.push_back( lsuma );
    phasea.push_back( isign ); // Temporary storage

    set<int> sorted( (*ideter).begin(), (*ideter).end() );
    copy (sorted.begin(),sorted.end(),(*ideter).begin());
  }
  
  // Beta orbitals next
  
  for( ideter=ioutfg.beta_spins.begin(); ideter != ioutfg.beta_spins.end(); ++ideter ) {
    
    int norbs = (*ideter).size();
    
    vector<int> buffer(nbf,0);
    int lsumb = 0; 
    int isign = 0;
    
    for(it = (*ideter).begin(); it != (*ideter).end(); ++it ) {
      int orb = (*it);
      lsumb += orb+1; 
      buffer[orb] = 1;
      ++num_processed;
      
      if ( orb != nbf-1 ) {
	for( int i=orb+1; i<nbf; ++i) {
	  isign += buffer[i];
	}
      }
      //  cout << "for B org " << orb << " FInal isign is " << isign << " mod 2 is " << isign % 2 << endl;
    }
    
    ioutfg.sumsb.push_back( lsumb );
    phaseb.push_back( isign ); // Temporary storage

    set<int> sorted( (*ideter).begin(), (*ideter).end() );
    copy (sorted.begin(),sorted.end(),(*ideter).begin());
  }
  
  // Assemble the phases
  
  if ( phasea.size() == phaseb.size() ) { 
    vector<int>::iterator ia,ib;
    
    ib=phaseb.begin(); 
    for(ia=phasea.begin(); ia != phasea.end(); ++ia){
      int sgna=(*ia);
      int sgnb=(*ib);
      ioutfg.phases.push_back( (sgna+sgnb) % 2 );
      ++ib;
    }
  } else {
    cerr << "Erroneous sizes: phase a <> phase b " << phasea.size() << " " << phaseb.size() << endl;
    GA::Terminate();
  }
  
#if 0
  cout << "Total number of phases computed is " << ioutfg.phases.size() << endl;
  cout << "Total number of sumsa computer is " << ioutfg.sumsa.size();
#endif
  
  ioutfg.completed_state = 1;

  return( num_processed );
}

// Collect the thee important assemble/distribution calls into one
/* These calls ensure that even if a node hasn't read or processed any 
   configuratyions/determinants that they all know what the basic job
   conditions are such as maxspatials, total sefs, etc.
*/ 
// Collective.

int PsociDeterminants::distributeParameters()
{
    assembleGlobalNumSpatials();
    assembleGlobalMaxDetSef();
    assembleGlobalJobParameters();
    return( 0 );
}

void PsociDeterminants::copyInternalMapToVector()
{
    local_configs->copyInternalMapToVector();
}




