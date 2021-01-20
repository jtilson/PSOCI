/**************************************************************************************

* Copyright (c) 2010,2011 RENCI.
* All rights reserved. This program and the accompanying materials
* MAY BE available under the terms of the RENCI Open Source License
* UNC at Chapel Hill which accompanies this distribution, and is available at
* http://www.renci.org/resources/open-source-software-license


 Classes: 

 Description: A series of functions to perform interpretation of the CI results

 History:

 

**************************************************************************************/
/**
 *   @file PsociCIAnalysis.C
 *
 */

#include <list>
#include "PsociCIAnalysis.hpp"

using namespace std;

/* We actually do not need the GA hamiltonian arrays cimat, icicol, and number (or supps) but we do need
   many of the parameters that charaterize the problem.The hamiltonian.destoy() is the proper method to use.
*/

// Currently only a single processor does the work
// All others are sent back to the calling program
/* Added a new piece. It may be possible that for some calculations (such as a VFCI) one may want to extract the list of configurations
   Most relevant to a select number of roots. In this way one could consider winnowing the reference space from which the CI is
   ultimately performed.  THe API for this method is not perfect  the number of roots to sumover AND the threshold are somewhat clunky

   leverage thresh since we sort of do this already
*/

int analyzeEnergyContributions( pair<double,double> & thresh, const double fcore,  const vector<int> & l_nsef, const vector<double> & evals, 
				const vector<double> & diags, GA::GlobalArray * g_v, PsociGAbasis & vector_set, 
				PsociGAhamiltonian & hamiltonian, PsociGADeterminants & deters)
{
  if (GA::nodeid() != 0 ) return(0);

  int retvalue = 1; // guilty

  const double dZero = 0.0;
  const double dOne = 1.0;

  cout << "/// Begin basis processing of the set of eigenvalues. "<<endl;
  cout << "    compute the energy contributions of the important spatial configurations" << endl;
  cout << "    to all eigenvectors: analyzeEnergyContributions: Roots not necessarily converged " << endl;
  cout << "    Check previously dumped rnorms for convergence check /// " << endl << endl;

  bool DumpRefs=false;
  if ( thresh.second > 0.0 ) { 
    cout << "RESER DumpRefs " << endl;
    DumpRefs=true;
  }

  set<string> refs; // in theory for use by CGDBGIN

  if ( DumpRefs ) {
  if ( thresh.second < thresh.first ) {
    cout << "thresh.second < thresh.first: resetting thresh.second to " << thresh.first << endl;
    thresh.second = thresh.first;
  }
  cout << "Will dump configurations in CGDBGIN format up through a threshhold of " << thresh.second << endl;
  }

  if ( l_nsef.size() != deters.fetchMaxSpatials() ) {
    GA::error("l_nsef != deters.fetchMaxSpatials()" ,-1);
  }
  
  // Grab vector globals
  
  int numroots =   vector_set.CurrentNumSoughtRoots();

// Grab the determinant globals
  int l_ksym = deters.fetchGlobalSymmetry();
  int l_maxSpatials = deters.fetchMaxSpatials();
  int l_maxsef = deters.fetchGlobalMaxSefs();
  int l_maxsefperconf = deters.fetchGlobalMaxSefPerSpatial();
  int l_maxdetperconf =  deters.fetchGlobalMaxDetPerSpatial();
  int l_nelec = deters.fetchNumElectrons();
  int l_nbfn = deters.fetchNumBasisFunctions();
  
  
  // Get qualities of the distributed vectors
  int type;
  int ld;
  int ndim, dims[2];
  
  g_v->inquire( &type, &ndim, dims );
  
  int lo[2], hi[2];
  
  int isefstart, isefstop;
  double vect2;
  JOUTFG fetch; 
  int index;
  
  string configs; // String to hold a representation of the spatial configuration
  double del;
  
  vector<double> vect( dims[0], dZero );
  double ebase = dZero;
  double energy_diff;
  
  const int wdth = l_nbfn+2; // Used in output formatting: configs

  for(int iroot=0; iroot < numroots; ++iroot ) 
    {
      if ( iroot==0 ) ebase = evals[iroot];
      
      //override data-parallel lo/hi for testing
      lo[0] = 0;
      hi[0] = dims[0]-1;
      lo[1] = iroot;
      hi[1] = iroot;
      ld = 1;
      g_v->get( lo, hi, &vect[0], &ld ); 
      
//Process sef related terms for the given ROOT iroot

// a) compute the energy contribution for spatial configurations of the ith root.

      vector<double> deltae( l_maxSpatials, dZero );
      for(int itotfg=0; itotfg< l_maxSpatials; ++itotfg ) 
	{
	  deltae[itotfg]=dZero;
	  isefstart = l_nsef[ itotfg ]; 
	  ( itotfg + 1 < l_maxSpatials ) ? isefstop = l_nsef[ itotfg + 1 ]: isefstop = l_maxsef;
	  for(int i=isefstart; i< isefstop; ++i ) 
	    {
	      vect2 = vect[i]*vect[i];
	      (vect2 == dOne) ? deltae[itotfg] = evals[iroot]: deltae[itotfg] -= vect2*(diags[i]-evals[iroot])/( dOne-vect2 );
	    }
	}

      // Now we want to sort the contributions by energy      
      pair<int,double> temppair;
      list<pair<int,double> > sortlist;
      for(int itotfg=0; itotfg< l_maxSpatials; ++itotfg )
	{
	  temppair.first=itotfg;
	  temppair.second=deltae[itotfg];
	  sortlist.push_back( temppair ); // comparison is based on abs() value
	}
 
      // Keep deltae array fopr the statistical methods

      sortlist.sort( compare_pairs ); //sort in reverse order based on abs(pair.second) value

      int count = 0;
      double value=dZero;
      
      //Truncated and sorted list:  Find first value < threshold, skip the rest ( they are ordered )
      list<pair<int,double> >::iterator it;
      for(it=sortlist.begin(); it!= sortlist.end(); ++it)
	{
	  value = (*it).second;
	  if ( abs(value) < thresh.first ) break;
	  ++count;
	}
      sortlist.resize( count );
      
       
/* Write out contributions for (abs) energy > threshold
   Exploit the fact that the sortlist is SORTED by energy 
   No need to go further down the list once below thresh.

   Need to fetch data for the determinants 
*/

      int l_nbfn = deters.fetchNumBasisFunctions();
      
      energy_diff = (evals[iroot] - ebase) * AUtoEV;
      cout << endl;
      cout << setfill(' ')<<setw(50)<<"Eigenstate"<< "\t\t"<< iroot+1 << endl;
      cout << setiosflags(ios::right);
      cout << setiosflags(ios::fixed);
      cout << setw(60)<<"Valence energy (au) "<<setw(12)<<setprecision(8)<<evals[iroot]<<endl;
      cout << setw(60)<<"Core energy (au)    "<<setw(12)<<setprecision(8)<<fcore<<endl;
      cout << setw(60)<<"Total energy (au)   "<<setw(12)<<setprecision(8)<<evals[iroot]+fcore<<endl;
      cout << setw(60)<<"E(i) - E(1) (eV)    "<<setw(12)<<setprecision(8)<<energy_diff<<endl;
      cout << endl;
      

      cout <<setw(9)<<"nsp"<<setw(1+wdth/2)<<"Spatial configuration"<<setw(wdth/2)<<" "<<setw(3)<<"sym"<<" "<<setw(9)<<"ndbg";
      cout <<" "<<setw(9)<< "vect[i]"<<" "<<setw(9)<<"dbg dE"<<"  "<<setw(9)<<"cnfg dE"<<endl<< endl;

#ifdef DIRECT
      vector<COMPRESS_SIZE> imap( l_nbfn+2  );
#endif

/* At this point we loop over all confs in sort. But thresh.second >= thresh.first so we
   can simply fetch refs here
*/

      for(it=sortlist.begin(); it!=sortlist.end(); ++it) { // new loop over configs
	
	index = (*it).first;

#ifdef DIRECT
      deters.fetchOrbMap( index+1, imap);
#endif

#ifdef DIRECT
        deters.computeConfigsGA( index+1, imap, fetch );
#else
	deters.fetchAndUnpackDeterminantData( index+1, fetch );
#endif
	deters.constructConfigurationData( fetch, configs ) ; // a conveninence for the user

// Currently inserts for ALL roots - may want to prune that in the future

      if ( DumpRefs && abs((*it).second) < thresh.second ) {
         refs.insert( configs );
      }

// b) compute the sef-sef contribution for particular spatial configuration for the ith root.

	isefstart = l_nsef[ index ]; 
	( index + 1 < l_maxSpatials ) ? isefstop = l_nsef[ index + 1 ]: isefstop = l_maxsef;
	for(int i=isefstart; i< isefstop; ++i ) 
	  {
	    vect2 = vect[i]*vect[i];
	    (vect2 == dOne) ? del = evals[iroot]: del = vect2*(evals[iroot] - diags[i])/( dOne-vect2 );
	    if ( i == isefstart ) {

	      cout <<setw(9)<<(*it).first+1<<setw(wdth)<<configs<<" "<<setw(3)<<fetch.jsym<<" "<<setw(9)<<i+1;
	      cout <<" "<<setprecision(6)<<setw(10)<< vect[i]<<" "<<setw(10)<<del<<" "<<setw(10)<<(*it).second<<endl;
	      
	    } else {
	      cout << setw(9+3+2+wdth)<<" "<<setw(9)<<i+1<<" "<<setprecision(6)<<setw(10)<<vect[i]<<" "<<setw(10)<<del<<endl;
	    }
	  }  
      }
      
// Statistical analysis and printout
// Return value of no use for now
      
      if ( energyDistributionSpatial( thresh.first, deltae ) != 0 ) {
	cout<<" Warning: Distribution method failed: skipping: root " << iroot << endl;
      }
      
    } // iroot

// Dump refs in CGDBGIN format if requested

  //string excitationlevel="100";
  //string excitationlevel="200";
  //cout << "at the dump step " << endl;
  //if ( DumpRefs ) {
     //outputSelectedRefstoCGDBG( refs, excitationlevel );
  //}

  refs.clear();
  retvalue = 0;
  return( retvalue ); 
}

// Compare based on abs(deltae) value
int energyDistributionSpatial( const double thresh, const vector<double> & deltae )
{
  double dZero = 0.0;
  int retvalue = 1; //guilty !

//   Default binning
  double bindata[]={0.00000001, 0.0000001, 0.000001, 0.00001, 0.0001, 0.001, 0.01, 0.1, 1.0, 10.0 };
  vector<double> bins(bindata, bindata+(sizeof(bindata)/sizeof(double)) );

//Ensure the data are sorted

  list<double> sortlist( deltae.begin(), deltae.end() );
  sortlist.sort( compare_values );
  
  int numConfigs = sortlist.size();
  
  double value, absvalue;
  
  int numCut=0;
  double sumCut = dZero;
  
  sumCut = dZero;
  numCut = 0;
  
  int binindex=0;
  list<double>::reverse_iterator rlit;

  cout << endl;
  cout << "\t\tCumulative contributions of spatial configurations to the current root " << endl;
  cout << "\t\t"<<setfill('-')<<setw(70)<<"-"<<endl;

/* Add a minor correction: what happens if a config has a contribution > 1?
   simply identify it and print out the number of them
*/
  for (rlit=sortlist.rbegin(); rlit != sortlist.rend(); ++rlit) //already sorted start from the back
    {
      value = (*rlit) ;
      absvalue = abs( value );
      
    checkloop:
      if ( absvalue >= bins[binindex] ) {

	cout << "\t\t"<<"Cutoff (au) "<<scientific<<setprecision(1)<<bins[binindex]<< '\t'<<" Num Spatials ";
	cout << setfill(' ')<<fixed<<setprecision(6)<< setiosflags(ios::left)<<setw(7)<<numCut<<" Energy sum of contributions (au) ";
	cout << setprecision(8)<< setiosflags(ios::left)<<sumCut<<endl;

	sumCut = dZero;
	numCut = 0;
        if ( binindex < bins.size() ) {
           ++binindex;
        } else {// Since the data are sorted once we hit the end  we can safely stop
          break;
        } 
	goto checkloop;
      }

      ++numCut;
      sumCut += value;
    }

  cout << "\t\t"<<"Cutoff (au) "<<scientific<<setprecision(1)<<bins[binindex]<< '\t'<<" Num Spatials ";
  cout << fixed<<setprecision(6)<< setiosflags(ios::left)<<setw(7)<<numCut<<" Energy sum of contributions (au) ";
  cout << setprecision(8)<< setiosflags(ios::left)<<sumCut<<endl;
  cout << endl<<endl;

  retvalue = 0;
  return( retvalue );
}

// Sort from abs(value) highest to lowest) - zeros at the bottom
bool compare_pairs( pair<int,double> left, pair<int,double> right )
{
  if (abs(left.second) > abs(right.second) ) return( true );
  return ( false );
}

bool compare_values( double left, double right )
{
  if (abs(left) > abs(right) ) return( true );
  return ( false );
}

/* A somewhat hacky way to retrieve a set of important spatials as determined
   by a lower level PSOCI calculation. E.g., print a set of configs that were globally important
   for a VFCI. Then perform +S or +S+D on that set.
  
   The plan here is to dump this set to a file for cut/past to a cgdbgin file
*/
void outputSelectedRefstoCGDBG( set<string> & refs, string & excitationlevel )
{
     if ( refs.size() == 0 ) {
        cout << "Empty refs set: Nothing to dump: return" << endl;
        return;
     }

     ofstream ofile;
     const string filename="cgdbgin.refs";
     ofile.open( filename.c_str() ); //rewrite any existing file
     if (!ofile.is_open()) {
       cout << "Failed to open refs dump file " << filename << endl;
       return;
     }
     set<string>::iterator sit;
// num basic configs , num configs specifically put on the list
       ofile <<" num basic configs , num configs specifically put on the list follows" << endl;
       ofile <<" these are entries 4 and 5, respectively on cgdbgin line 2 (first line after title) " << endl;
       ofile <<" note also the following list INCLUDEs the trailing 0 required by cgdbgin " << endl;
       ofile << endl;
       ofile << refs.size() << " "<< refs.size() << endl;
       ofile << endl;

//This places refs on the cgdbgin list
     for( sit=refs.begin(); sit != refs.end(); ++sit ) {
        ofile << (*sit) << endl;
     }
//This places refs with associated excitation level on the list

     for( sit=refs.begin(); sit != refs.end(); ++sit ) {
        ofile << (*sit) << excitationlevel << endl;
     }

// Place the trainling '0' that cgdbgin requires
     ofile << '0' << endl;

     refs.clear();
     ofile.close();
     cout << "Dump of selected references to CGDBGIN completed " << endl;
     return;
}
 






