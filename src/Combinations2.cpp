

/**
 * @file Combinations.cpp
 * Find all the combinations of n things taken k at a time.
 */
 

#include "Combinations2.hpp"

namespace gpstk {


// ----------------------------------------------------------------------------
// Combinations.cpp
// Find all the Combinations of n things taken k at a time.


/// Class Combinations will compute C(n,k), all the Combinations of n things
/// taken k at a time (where k <= n).
/// Let n 'things' be indexed by i (i=0...n-1), e.g. stored in an array of length n.
/// This class computes C(n,k) as sets of k indexes into the 'things' array.
/// These indexes are accessible via member functions Selection() or isSelected().
/// Next() computes the next combination until there are no more (when it returns -1).

// ----------------------------------------------------------------------------

      /// Default constructor
   Combinations::Combinations(void)
      throw()
   {
      nc = n = k = 0;
      Index = NULL;
   }
   
      /// Constructor for C(n,k) = Combinations of n things taken k at a time (k <= n)
      /// @throw on invalid input (k>n).
   Combinations::Combinations(int N, int K)
      throw(gpstk::Exception)
   {
      try
      {
         init(N,K);
      }
      catch(gpstk::Exception& e)
      {
         GPSTK_RETHROW(e);
      }
   }
   
      /// copy constructor
   Combinations::Combinations(const Combinations& right)
      throw()
   {
      *this = right;
   }
   
      /// destructor
   Combinations::~Combinations(void)
   {
      if (Index)
         delete[] Index;
      Index = NULL;
   }
   
      /// Assignment operator.
   Combinations& Combinations::operator=(const Combinations& right)
      throw()
   {
      this->~Combinations();
      init(right.n,right.k);
      nc = right.nc;
      for (int j=0; j<k; j++)
      {
         Index[j] = right.Index[j];
      }
      return *this;
   }
   
      /// Return index i (0 <= i < n) of jth selection (0 <= j < k);
      /// if j is out of range, return -1.
   int Combinations::Selection(int j)
      throw()
   {
      if (j < 0 || j >= k)
         return -1;
      return Index[j];
   }
   
      /// Return true if the given index j (0 <= j < n) is
      /// currently selected (i.e. if j = Selection(i) for some i)
   bool Combinations::isSelected(int j)
      throw()
   {
      for (int i=0; i<k; i++)
      {
         if (Index[i] == j)
            return true;
      }
      return false;
   }


   
   

// ----------------------------------------------------------------------------

      /// Compute the next combination, returning the number of Combinations computed
      /// so far; if there are no more Combinations, return -1.
int Combinations::Next(void)
   throw()
{
   if (k < 1)
      return -1;
   if (Increment(k-1) != -1)
      return ++nc;
   return -1;
}

      /// Recursive function to increment Index[j].
int Combinations::Increment(int j)
   throw()
{
   if (Index[j] < n-k+j) {        // can this index be incremented?
      Index[j]++;
      for (int m=j+1; m<k; m++)
      {
         Index[m]=Index[m-1]+1;
      }
      return 0;
   }
      // is this the last index?
   if (j-1 < 0) return -1;
      // increment the next lower index
   return Increment(j-1);
}

      /// The initialization routine used by constructors.
      /// @throw on invalid input (k>n or either n or k < 0).
void Combinations::init(int N, int K)
   throw(gpstk::Exception)
{
   if (K > N || N < 0 || K < 0)
   {
      gpstk::Exception e("Combinations(n,k) must have k <= n, with n,k >= 0");
      GPSTK_THROW(e);
   }
   
   if (K > 0)
   {
      Index = new int[K];
      if (!Index)
      {
         gpstk::Exception e("Could not allocate");
         GPSTK_THROW(e);
      }
   }
   else
      Index = NULL;
   
   nc = 0;
   k = K;
   n = N;
   for (int j=0; j<k; j++)
   {
      Index[j]=j;
   }
}
} // end namespace combinations
