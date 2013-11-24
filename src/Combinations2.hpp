
#pragma ident "$Id$"

/**
 * @file Combinations.hpp
 * Find all the combinations of n things taken k at a time.
 */
 
#ifndef COMBINATIONS_INCLUDE
#define COMBINATIONS_INCLUDE



// ----------------------------------------------------------------------------
// Combinations.hpp
// Find all the Combinations of n things taken k at a time.

#include <vector>
#include <ostream>
#include "GNSSconstants.hpp"
#include "CommonTime.hpp"
#include "SatID.hpp"
#include "Matrix.hpp"
#include "RinexObsHeader.hpp"
#include "XvtStore.hpp"
#include "TropModel.hpp"
#include "Combinations2.hpp"
#include "Exception.hpp"

#include "MathBase.hpp"
#include "EphemerisRange.hpp"
#include "GPSEllipsoid.hpp"
#include "GPSWeekSecond.hpp"


namespace gpstk {

/// Class Combinations will compute C(n,k), all the Combinations of n things
/// taken k at a time (where k <= n).
/// Let n 'things' be indexed by i (i=0...n-1), e.g. stored in an array of length n.
/// This class computes C(n,k) as sets of k indexes into the 'things' array.
/// These indexes are accessible via member functions Selection() or isSelected().
/// Next() computes the next combination until there are no more (when it returns -1).

// ----------------------------------------------------------------------------
class Combinations
{

public:

      /// Default constructor
   Combinations(void)
      throw();
   
      /// Constructor for C(n,k) = Combinations of n things taken k at a time (k <= n)
      /// @throw on invalid input (k>n).
   Combinations(int N, int K)
      throw(gpstk::Exception);
   
      /// copy constructor
   Combinations(const Combinations& right)
      throw();
   
      /// destructor
   ~Combinations(void);
   
      /// Assignment operator.
   Combinations& operator=(const Combinations& right)
      throw();
   
      /// Compute the next combination, returning the number of Combinations computed
      /// so far; if there are no more Combinations, return -1.
   int Next(void) throw();
   
      /// Return index i (0 <= i < n) of jth selection (0 <= j < k);
      /// if j is out of range, return -1.
   int Selection(int j)
      throw();
   
      /// Return true if the given index j (0 <= j < n) is
      /// currently selected (i.e. if j = Selection(i) for some i)
   bool isSelected(int j)
      throw();

private:

      /// The initialization routine used by constructors.
      /// @throw on invalid input (k>n or either n or k < 0).
   void init(int N, int K)
      throw(gpstk::Exception);
   
      /// Recursive function to increment Index[j].
   int Increment(int j) throw();
   
   int nc;         ///< number of Combinations computed so far
   int k,n;        ///< Combinations of n things taken k at a time, given at c'tor
   int* Index;     ///< Index[j] = index of jth selection (j=0...k-1; I[j]=0...n-1)

};

} // end namespace combinations
#endif // COMBINATIONS_INCLUDE
