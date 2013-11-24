#pragma ident "$Id$"

/**
 * @file PRSolution2.hpp
 * Autonomous pseudorange navigation solution, including RAIM2 algorithm
 */
 
#ifndef PRS_POSITION_SOLUTION_LD_HPP
#define PRS_POSITION_SOLUTION_LD_HPP

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
//  Copyright 2004, The University of Texas at Austin
//
//============================================================================

#include <vector>
#include <ostream>
#include "GNSSconstants.hpp"
#include "CommonTime.hpp"
#include "SatID.hpp"
#include "Matrix.hpp"
#include "RinexObsHeader.hpp"
#include "XvtStore.hpp"
#include "TropModel.hpp"
#include "PRSolution2.hpp"
#include "Combinations2.hpp"

namespace gpstk
{
   /** @defgroup GPSsolutions GPS solution algorithms and Tropospheric models */
   //@{
 
   /** This class defines an interface to routines which compute a position
    * and time solution from pseudorange data, with a data editing algorithm
    * based on Receiver Autonomous Integrity Monitoring (RAIM) concepts.
    * RAIM ref. "A Baseline GPS RAIM Scheme and a Note on the Equivalence of
    * Three RAIM Methods," by R. Grover Brown, Journal of the Institute of
    * Navigation, Vol. 39, No. 3, Fall 1992, pg 301.
    */
   class PRSolution3 : public PRSolution2
   {
   public:
   /*
         /// Constructor
      PRSolution3() throw() :
         RMSLimit(6.5), SlopeLimit(1000.), Algebraic(false),
         ResidualCriterion(true), ReturnAtOnce(false), NSatsReject(-1),
         Debug(false), pDebugStream(&std::cout), MaxNIterations(10),
         ConvergenceLimit(3.e-7), Valid(false)
      {};
*/
      /** Compute a position/time solution, given satellite PRNs and pseudoranges
       *  using a RAIM algorithm.
       * @param Tr          Measured time of reception of the data.
       * @param Satellite   std::vector<SatID> of satellites; on successful
       *                    return, satellites that were excluded by the algorithm
       *                    are marked by a negative 'prn' member.
       * @param Pseudorange std::vector<double> of raw pseudoranges (parallel to
       *                    Satellite), in meters.
       * @param Eph         gpstk::EphemerisStore to be used in the algorithm.
       *
       * @return Return values:
       *  2  solution is found, but it is not good (RMS residual exceed limits)
       *  1  solution is found, but it is suspect (slope is large)
       *  0  ok
       * -1  algorithm failed to converge
       * -2  singular problem, no solution is possible
       * -3  not enough good data (> 4) to form a (RAIM) solution
       *     (the 4 satellite solution might be returned - check isValid())
       * -4  ephemeris is not found for one or more satellites
       */
      int RAIM2Compute(const CommonTime& Tr,
                      std::vector<SatID>& Satellite,
                      const std::vector<double>& Pseudorange,
                      const XvtStore<SatID>& Eph,
                      TropModel *pTropModel)
         throw(Exception);

         /// Return the status of solution
      bool isValid()
         const throw() { return Valid; }

      // output:

      // The bias estimate of every SV pseudorange at this epoch
      //Vector<double> BiasEstimates;

      // for every combination of 4 satellites, store the solution in a solution matrix
      std::vector< Vector<double> > ComboSolution;
      std::vector< std::vector<bool> > ComboSolSat; 

   }; // end class PRSolution3

   //@}

} // namespace gpstk

#endif
