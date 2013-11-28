#pragma ident "$Id$"

//============================================================================
//
//  This file is part of GPSTk, the GPS Toolkit.
//
//  The GPSTk is free software; you can redistribute it and/or modify
//  it under the terms of the GNU Lesser General Public License as published
//  by the Free Software Foundation; either version 2.1 of the License, or
//  any later version.
//
//  The GPSTk is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU Lesser General Public License for more details.
//
//  You should have received a copy of the GNU Lesser General Public
//  License along with GPSTk; if not, write to the Free Software Foundation,
//  Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110, USA
//
//  Copyright 2009, The University of Texas at Austin
//
//============================================================================


// GPSTk example program #4


   // First, let's include Standard Template Library classes
#include <string>
#include <vector>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <fstream>

   // Classes for handling observations RINEX files (data)
#include "Rinex3ObsHeader.hpp"
#include "Rinex3ObsData.hpp"
#include "Rinex3ObsStream.hpp"

   // Classes for handling satellite navigation parameters RINEX
   // files (ephemerides)
#include "Rinex3NavHeader.hpp"
#include "Rinex3NavData.hpp"
#include "Rinex3NavStream.hpp"

   // Classes for handling RINEX files with meteorological parameters
#include "RinexMetBase.hpp"
#include "RinexMetData.hpp"
#include "RinexMetHeader.hpp"
#include "RinexMetStream.hpp"

   // Class for handling tropospheric models
#include "TropModel.hpp"

   // Class for storing "broadcast-type" ephemerides
#include "GPSEphemerisStore.hpp"

   // Class for handling RAIM
#include "PRSolution3.hpp"

   // Class defining GPS system constants
#include "GNSSconstants.hpp"

//vplot
#include "EPSImage.hpp"
#include "SurfacePlot.hpp"
#include "Frame.hpp"
#include "BorderLayout.hpp"
#include "HLayout.hpp"
#include "VLayout.hpp"


using namespace std;
using namespace gpstk;

using namespace vdraw;
using namespace vplot;


long power2(int exponent) {
  return ((long)1 << exponent);
}

void debugPrintCombo(vector< bool > sats, vector < SatID > sids) {
               // now print all the little solutions
               // debug out FIXME print solutions of all combos
                 cerr << "Satellites ";
		 for (int g=0; g<sats.size(); g++)
		      cerr << sats[g] << " ";
		 cerr << endl;
		 for (int g=0; g<sids.size(); g++)
		      cerr << sids[g].id << " ";
		 cerr << endl;
}

void debugPrintVec(const char * label, Vector<double> vec) {
                 cerr << label << " " 
		      << vec[0] << " "
		      << vec[1] << " "
		      << vec[2] << " "
		      << vec[3] << " "
		      << endl;
		      }

class ComboContainer {
  public:
    //ComboContainer();
    
    //vector< map< SatID, Vector<double> > > prnVec;
    
    // List of all satellite prns
    vector< SatID > prns; // was SatID
    vector < int > sol_ids;
    // by epoch time series of combo solutions.
    vector< map< int, Vector< double > > > ComboSolutions;

    int pushSolID(int id) {
      std::vector< int >::iterator it;
      if ((it = std::find(sol_ids.begin(), sol_ids.end(), id))==sol_ids.end())
      {
        sol_ids.push_back(id);
        return 1;
      } else return 0;
    }

    void sortSolIDs() {
      std::sort(sol_ids.begin(), sol_ids.end());
    }

    // from the above list of prns we determine all combination IDs
    int getComboID(vector< SatID > in_prns, vector< bool > use_sv) {
      // are all prns in the prns vector?
      long id=0;
      for (int i=0;i<in_prns.size();i++) {
        std::vector< SatID >::iterator it;
	SatID abssat = in_prns[i];
	abssat.id = abs(abssat.id);
	int power;
        if ((it = std::find(prns.begin(), prns.end(), abssat))!=prns.end()){
	   power = distance(prns.begin(),it);
	} else {
	   power = prns.size();
	   prns.push_back(abssat);
	}
        if (use_sv[i]) id += power2(power);
      }
      pushSolID(id); // who cares about return value right now.
      return id;
    }

    int newEpoch() {
      ComboSolutions.push_back(map< int, Vector< double > >());
    }

    void storeSolution(int id, Vector< double > sol) {
      ComboSolutions[ComboSolutions.size()-1][id]=sol;
    }

    int getTimeSeries(int id, int epoch, int length, vector< Vector<double> > & ts) {
      ts.clear();
      // if an unbroken time series of this SV id exists, return 1.
      // otherwise return 0;

      if (ComboSolutions.size() < epoch || (epoch - length + 1) < 0) return 0;
      while (length > 0) {
        // get the solution of SV[id] for epoch
	if (ComboSolutions[epoch][id].size() == 0) return 0; // this creates an entry (thus memory management) but is faster than catching an exception. Blame STL for this mess.
	ts.push_back(ComboSolutions[epoch][id]);
	epoch--;
        length--;
      }
      return 1;
    }

    Vector< double > calculateDelta(vector< Vector< double > > vi, vector< Vector< double > > vj) {
      Vector<double> diffsum(vi[0]);
      for (int i=0; i<diffsum.size(); i++) diffsum[i]=0;
      for (int t=0;t<vi.size();t++){
        diffsum += (vi[t]-vj[t]);
      }
      diffsum = diffsum / (double)(vi.size());
      return diffsum;
    }

    // TODO make this work on positioning solution, not just one component.
    // for now, select component using index cp = {0,1,2 or 3};
    double calculateBeta(int cp, vector< Vector< double > > vi, double biasi, vector< Vector< double > > vj, double biasj) {
      int T = vi.size();
      double beta=0;
      for (int t=0;t<T;t++) {
        double term = (vi[t][cp] - biasi) - (vj[t][cp] - biasj);
	beta += term*term;
      }
      return (beta / (T-1));
    }


    // class variables
    vector< vector< Vector< double > > > all_ts;
    vector< int > use_sol_ids;

    Vector< double > bias;
    Vector< double > variance;

    int calculateBiases() {
      sortSolIDs();
      int T = 50;
      int unusedSols = 0;
      // for now just use the fourth component (receiver clock error)
      int cp=0;

      cerr << "Number of solution IDs = " << sol_ids.size() <<endl;

      
      // Now we can calculate the biases of a time series of length say 5
      for (int i=0;i<sol_ids.size();i++) {
        int id=sol_ids[i];
        vector< Vector<double> > ts;
        if (getTimeSeries(id, T, T, ts)) {
	  all_ts.push_back(ts);
	  use_sol_ids.push_back(id);
	} else {
	  unusedSols++;
	}
      }

      cerr << "Unused Solutions = " << unusedSols << endl;

      int dim=use_sol_ids.size();

      Matrix< Vector< double > > delta ( dim, dim);

      for (int i=0; i< dim; i++) {
        int id_i = use_sol_ids[i];
        for (int j=0; j< dim; j++) {
          int id_j = use_sol_ids[j];
	  // add this to the delta matrix;
          delta[i][j]=calculateDelta(all_ts[i], all_ts[j]);
	  //debugPrintVec("Delta " , delta[i][j]);
	}
      }

      // Deltas have been computed. Now it is time to solve the system of equations (3.3 in AI's notes)
      Matrix< double > A (dim + 1,dim + 1); // + 1 for lambda parameter
      Vector< double > b (dim+1);

      cerr << "Generate A" << endl;
      for (int k=0; k < dim; k++) {
        //cerr << "k=" << k << endl;
        A[k][k]=0;
        for (int i=0; i< dim; i++) {
	  if (i!=k) {
            A[k][i]=2.0/(delta[k][i][cp] * delta[k][i][cp]); // 2/(d(i,k)^2)
	    A[k][k]+=A[k][i]; // sum(2/(d(i,k)^2))
	  }
	}
	//A[k][k]=gpstk::sum<double>(A[k]); // sum(2/(d(i,k)^2))
	A[k][dim]=1; // \lambda
      }
      for (int i=0; i< dim; i++) A[dim][i]=1; // sum(b_i) == 0
      A[dim][dim]=0;

      cerr << "Generate b" << endl;
      // b vector
      for (int k=0; k < dim; k++) {
        b[k]=0;
        for (int i=0; i< dim; i++) {
          if (i!=k) b[k] += ((i >= k) ? 1.0 : -1.0) / delta[k][i][cp];
	}
	b[k] *= 2;
      }
      b[dim] = 0; // sum of b_i == 0

      cerr << "Solve bias=A^(-1) b" << endl;
      // now solve
      bias (dim+1);
      bias = gpstk::inverse(A)*b; // heaps slow.

      // can I print bias?
      cerr << "Bias solution " << endl;
      for (int k=0; k < dim+1; k++) cerr << bias[k] << endl;

      // now calculate the variances
      cerr << "Calculating variances" << endl;

      // TODO FIXME make this Matrix< Vector< double > >
      Matrix< double > Beta ( use_sol_ids.size(), use_sol_ids.size());

      for (int i=0; i< use_sol_ids.size(); i++) {
        int id_i = use_sol_ids[i];
        for (int j=0; j< use_sol_ids.size(); j++) {
          int id_j = use_sol_ids[j];
	  // Beta matrix entry
          Beta[i][j]=calculateBeta(cp, all_ts[i], bias[i], all_ts[j], bias[j]);
	  //debugPrintVec("Delta " , delta[i][j]);
	}
      }
      
      // now formulate a second linear system using equations in 3.7 in AI's notes

      cerr << "Generate A" << endl;
      for (int k=0; k < dim; k++) {
        A[k][k]=0;
        for (int i=0; i< dim; i++) {
	  if (i!=k) {
            A[k][i]=1.0/(Beta[i][k] * Beta[i][k]); // 2/(d(i,k)^2)
	    A[k][k]+=A[k][i]; // sum(2/(d(i,k)^2))
	  }
	}
	//A[k][k]=gpstk::sum<double>(A[k]); // sum(2/(d(i,k)^2))
	A[k][dim]=0.5; // \lambda
      }
      for (int i=0; i< dim; i++) A[dim][i]=1; // sum(v_i) == the complicated rhs
      A[dim][dim]=0;

      cerr << "Generate b" << endl;
      // b vector
      for (int k=0; k < dim; k++) {
        b[k]=0;
        for (int i=0; i< dim; i++) {
          if (i!=k) b[k] += 1.0 / Beta[i][k]; // already [cp]
	}
      }
      cerr << "Calculate expected values using biased solutions" << endl;
      // FIXME TODO make this a vector of expected position solutions
      vector< double > avg_sol(T,0);
      for (int t=0; t< T; t++) {
        //avg_sol.push_back(0); // new epoch
        for (int i=0; i< dim; i++) {
          avg_sol[t] += all_ts[i][t][cp];
	}
	avg_sol[t] /= dim;
      }
      cerr << "Final equation in 3.7" <<endl;
      // b_{N+1}= rhs of bottom of eqn 3.7 
      for (int i=0; i < dim; i++) {
        double sum_sq_dx = 0;
        for (int t=0; t< T; t++) {
	  double dx= all_ts[i][t][cp]-avg_sol[t];
	  sum_sq_dx += dx * dx;
	}
	b[dim] += sum_sq_dx;
      }
      b[dim] *= ((double)dim)/(T*(dim-1));



      cerr << "Solve variance=A^(-1) b" << endl;
      // now solve
      variance (dim+1);
      variance = gpstk::inverse(A)*b; // Unoptimised. Should use L*L decomposition.
      // can I print?
      cerr << "Variance solution " << endl;
      for (int k=0; k < dim+1; k++) cerr << variance[k] << endl;

      

      // how do we relate position solution variances to prn variances?

      // it seems that 11C3 = 165 which means a lot of combinations share the error of just one pseudorange
      // therefore back-propagating the bias of this solution to the bias on the pseduorange

      //cerr << "B " << endl;
      //for (int k=0; k < dim+1; k++) cerr << b[k] << endl;

      //cerr << "A[0] " << endl;
      //for (int k=0; k < dim+1; k++) cerr << A[k][0] << endl;


    }

    void bestSolution(ostream * outfile) {
      *outfile << "test" << endl;

      double reg = 0.5; // regularisation factor \lambda

      int dim = variance.size()-1; // lambda is last term
      int T = all_ts[0].size();

      vector< double > weight(dim);


      double totalvariance = 0;
      // weight = 1/variance/(1/total variance)
      for (int k=0;k<dim;k++) {
        totalvariance+= 1.0/(variance[k] + bias[k]*bias[k] + reg);
      }

      for (int k=0;k<dim;k++) {
        weight[k] = 1.0/(variance[k]+bias[k]*bias[k] + reg) / totalvariance;
      }

      // use these to weight entire solution (not just component we're estimating)
      Vector< Vector< double > > ws(T,Vector< double >(4, 0));
      for (int t=0;t<T;t++) {
        for (int k=0;k<dim;k++) {
	  for (int i=0;i<4;i++ ) {
	    ws[t][i]+= weight[k]*all_ts[k][t][i];
	  }
	}
      }

      // write it out.
      for (int t=0;t<T;t++) {
	for (int i=0;i<4;i++) {
	  *outfile << setprecision(15) << setw(20) << ws[t][i];
	}
	*outfile << endl;
      }


        
    }
};


void debugPrintCombos(PRSolution3& raimSolver) {
               // now print all the little solutions
               // debug out FIXME print solutions of all combos
               for (int q=0; q<raimSolver.ComboSolution.size(); q++) {
                 cerr << "Solution   " 
		      << raimSolver.ComboSolution[q][0] << " "
		      << raimSolver.ComboSolution[q][1] << " "
		      << raimSolver.ComboSolution[q][2] << " "
		      << raimSolver.ComboSolution[q][3] << " "
		      << endl;
                 cerr << "Satellites ";
		 for (int g=0; g<raimSolver.ComboSolSat[q].size(); g++)
		      cerr << raimSolver.ComboSolSat[q][g] << " ";
		 cerr << endl;
             }
}

int main(int argc, char *argv[])
{

      // Declaration of objects for storing ephemerides and handling RAIM
   GPSEphemerisStore bcestore;
   PRSolution3 raimSolver;

      // Object for void-type tropospheric model (in case no meteorological
      // RINEX is available)
   ZeroTropModel noTropModel;

      // Object for GG-type tropospheric model (Goad and Goodman, 1974)
      // Default constructor => default values for model
   GGTropModel ggTropModel;

      // Pointer to one of the two available tropospheric models. It points
      // to the void model by default
   TropModel *tropModelPtr=&noTropModel;

      // This verifies the ammount of command-line parameters given and
      // prints a help message, if necessary
   if( (argc < 3) || (argc > 4) )
   {
      cerr <<  "Usage:" << endl;
      cerr << "   " << argv[0]
           << " <RINEX Obs file>  <RINEX Nav file>  [<RINEX Met file>]"
           << endl;

      exit (-1);
    }

      // Let's compute an useful constant (also found in "GNSSconstants.hpp")
   const double gamma = (L1_FREQ_GPS/L2_FREQ_GPS)*(L1_FREQ_GPS/L2_FREQ_GPS);

   try
   {
      ComboContainer combolist;

         // Read nav file and store unique list of ephemerides
      Rinex3NavStream rnffs(argv[2]);    // Open ephemerides data file
      Rinex3NavData rne;
      Rinex3NavHeader hdr;

         // Let's read the header (may be skipped)
      rnffs >> hdr;

         // Storing the ephemeris in "bcstore"
      while (rnffs >> rne) bcestore.addEphemeris(rne);

         // Setting the criteria for looking up ephemeris
      bcestore.SearchNear();

         // If provided, open and store met file into a linked list.
      list<RinexMetData> rml;
        
      if( argc == 4 )
      {
         cerr << "Loading meteorological file" << endl;

         RinexMetStream rms(argv[3]);    // Open meteorological data file
         RinexMetHeader rmh;
            
            // Let's read the header (may be skipped)
         rms >> rmh;

         RinexMetData rmd;
            
            // If meteorological data is provided, let's change pointer to
            // a GG-model object
         tropModelPtr=&ggTropModel;

            // All data is read into "rml", a meteorological data
            // linked list
         while (rms >> rmd) rml.push_back(rmd);

      }  // End of 'if( argc == 4 )'

         // Open and read the observation file one epoch at a time.
         // For each epoch, compute and print a position solution
      Rinex3ObsStream roffs(argv[1]);    // Open observations data file

         // In order to throw exceptions, it is necessary to set the failbit
      roffs.exceptions(ios::failbit);

      Rinex3ObsHeader roh;
      Rinex3ObsData rod;

         // Let's read the header
      roffs >> roh;

         // The following lines fetch the corresponding indexes for some
         // observation types we are interested in. Given that old-style
         // observation types are used, GPS is assumed.
      int indexC1;
      try
      {
         indexC1 = roh.getObsIndex( "C1" );
      }
      catch(...)
      {
         cerr << "The observation file doesn't have C1 pseudoranges." << endl;
         exit(1);
      }

      int indexC2;
      try
      {
         indexC2 = roh.getObsIndex( "C2" );
      }
      catch(...)
      {
         indexC2 = -1;
      }

         // Defining iterator "mi" for meteorological data linked list
         // "rml", and set it to the beginning
      list<RinexMetData>::iterator mi=rml.begin();

         // Let's process all lines of observation data, one by one
      while( roffs >> rod )
      {

            // Find a weather point. Only if a meteorological RINEX file
            // was provided, the meteorological data linked list "rml" is
            // neither empty or at its end, and the time of meteorological
            // records are below observation data epoch.
         while( ( argc==4 )       &&
                ( !rml.empty() )  &&
                ( mi!=rml.end() ) &&
                ( (*mi).time < rod.time ) )
         {

            mi++;    // Read next item in list

               // Feed GG tropospheric model object with meteorological
               // parameters. Take into account, however, that setWeather
               // is not accumulative, i.e., only the last fed set of
               // data will be used for computation
            ggTropModel.setWeather( (*mi).data[RinexMetHeader::TD],
                                    (*mi).data[RinexMetHeader::PR],
                                    (*mi).data[RinexMetHeader::HR] );

         }  // End of 'while( ( argc==4 ) && ...'


         // Apply editing criteria
         if( rod.epochFlag == 0 || rod.epochFlag == 1 )  // Begin usable data
         {

            vector<SatID> prnVec;
            vector<double> rangeVec;

               // Define the "it" iterator to visit the observations PRN map. 
               // Rinex3ObsData::DataMap is a map from RinexSatID to
               // vector<RinexDatum>:
               //      std::map<RinexSatID, vector<RinexDatum> >
            Rinex3ObsData::DataMap::const_iterator it;

               // This part gets the PRN numbers and ionosphere-corrected
               // pseudoranges for the current epoch. They are correspondly fed
               // into "prnVec" and "rangeVec"; "obs" is a public attribute of
               // Rinex3ObsData to get the map of observations
            for( it = rod.obs.begin(); it!= rod.obs.end(); it++ )
            {

                  // The RINEX file may have C1 observations, but the current
                  // satellite may not have them.
               double C1( 0.0 );
               try
               {
                  C1 = rod.getObs( (*it).first, indexC1 ).data;
               }
               catch(...)
               {
                     // Ignore this satellite if C1 is not found
		     cerr << "No C1 Obs found for epoch" << endl;
                  continue;
               }

               double ionocorr( 0.0 );

                  // If there are C2 observations, let's try to apply the
                  // ionospheric corrections
               if( indexC2 >= 0 )
               {

                     // The RINEX file may have C2 observations, but the
                     // current satellite may not have them.
                  double C2( 0.0 );
                  try
                  {
                     C2 = rod.getObs( (*it).first, indexC2 ).data;
                  }
                  catch(...)
                  {
                        // Ignore this satellite if C1 is not found
                     continue;
                  }

                     // Vector 'vecData' contains RinexDatum, whose public
                     // attribute "data" indeed holds the actual data point
                  ionocorr = 1.0 / (1.0 - gamma) * ( C1 - C2 );

               }

                  // Now, we include the current PRN number in the first part
                  // of "it" iterator into the vector holding the satellites.
                  // All satellites in view at this epoch that have C1 or C1+C2
                  // observations will be included.
               prnVec.push_back( (*it).first );

                  // The same is done for the vector of doubles holding the
                  // corrected ranges
               rangeVec.push_back( C1 - ionocorr );

                     // WARNING: Please note that so far no further correction
                     // is done on data: Relativistic effects, tropospheric
                     // correction, instrumental delays, etc.

            }  // End of 'for( it = rod.obs.begin(); it!= rod.obs.end(); ...'

               // The default constructor for PRSolution2 objects (like
               // "raimSolver") is to set a RMSLimit of 6.5. We change that
               // here. With this value of 3e6 the solution will have a lot
               // more dispersion.
            raimSolver.RMSLimit = 3e6;

               // In order to compute positions we need the current time, the
               // vector of visible satellites, the vector of corresponding
               // ranges, the object containing satellite ephemerides, and a
               // pointer to the tropospheric model to be applied
            raimSolver.RAIM2Compute( rod.time,
                                    prnVec,
                                    rangeVec,
                                    bcestore,
                                    tropModelPtr );

               // Note: Given that the default constructor sets public
               // attribute "Algebraic" to FALSE, a linearized least squares
               // algorithm will be used to get the solutions.
               // Also, the default constructor sets ResidualCriterion to true,
               // so the rejection criterion is based on RMS residual of fit,
               // instead of RMS distance from an a priori position.

               // If we got a valid solution, let's print it

            if( raimSolver.isValid() )
            {
                  // Vector "Solution" holds the coordinates, expressed in
                  // meters in an Earth Centered, Earth Fixed (ECEF) reference
                  // frame. The order is x, y, z  (as all ECEF objects)
               cout << setprecision(12) << raimSolver.Solution[0] << " " ;
               cout << raimSolver.Solution[1] << " " ;
               cout << raimSolver.Solution[2] << " " ;
               cout << raimSolver.Solution[3];

               cout << endl ;

               //debugPrintCombos(raimSolver); 


               // now push these results into
	       // some structure to hold a time series for each
	       // satellite prn

               combolist.newEpoch();

               // get the ID of this combo and store the solution
	       for (int g=0; g<raimSolver.ComboSolSat.size(); g++)
	          combolist.storeSolution(combolist.getComboID(prnVec,raimSolver.ComboSolSat[g]),raimSolver.ComboSolution[g]);


            }  // End of 'if( raimSolver.isValid() )'
	    else {
	      cout << "Invalid solution " << endl;
	    }

         } // End of 'if( rod.epochFlag == 0 || rod.epochFlag == 1 )'

      }  // End of 'while( roffs >> rod )'
      cerr << combolist.ComboSolutions.size() << " combination solutions stored" << endl;

      combolist.calculateBiases();
      
      ostream * outfile;
      string filename("bestoutput");
      outfile = new fstream((filename+".txt").c_str(), fstream::out);
      combolist.bestSolution(outfile);
      

   }
   catch(Exception& e)
   {
      cerr << e << endl;
   }
   catch (...)
   {
      cerr << "Caught an unexpected exception." << endl;
   }

   exit(0);

}  // End of 'main()'

