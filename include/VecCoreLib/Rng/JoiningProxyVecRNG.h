#ifndef SIMPLE_LOOP_PROXY_VECRNG_H
#define SIMPLE_LOOP_PROXY_VECRNG_H

#include "RngDefs.h"
#include "VecRNG.h"

namespace vecRng {
inline namespace VECRNG_IMPL_NAMESPACE {


template <typename BaseVecRNGType,
          typename BaseScalarRNGType,
          typename BackendT>
struct RNG_traits<JoiningProxyVecRNG<BaseVecRNGType, BaseScalarRNGType, BackendT>>
{
  // fCg[row=MRG::vsize][column=VectorSize<Double_v>()]
  struct StateJoined { typename BaseVecRNGType fVectorRNGs ; }
  using State_t = StateJoined;
};
   
template <typename BaseVecRNGType,
          typename BaseScalarRNGType,
          typename BackendT>
   class JoiningProxyVecRNG :
     // public BaseVecRNGType /* <BackendT> */    // Type 1: derive & extend
     public VecRNG< JoiningProxyVecRNG<BaseVecRNGType, BaseScalarRNGType, BackendT> >  // Type 2: Use 'has'
{
public:
   // using BaseRNGState_t= typename RNG_traits<BaseVecRNGType> /*<BackendT>>*/ ::State_t;
   using BaseRNGState_scalar= typename RNG_traits<BaseVecRNGType> /*<ScalarBackend>>*/ ::State_t;
   using BaseRNGState_vector= typename RNG_traits<BaseScalarRNGType> /*<VectorBackend>>*/ ::State_t;
   
   /** @brief Constructor: 'create' vector instance from a set of sequential states */
   JoiningProxyVecRNG( BaseState_scalar& trackPrng[/*VecLen*/], int /*numVariatesExpected*/ );

   /** @brief Split apart this class - ensure that the state of its parts is fed back */
   //          - i.e. the state gets copied back to each 'track' sequential RNG state
   ~JoiningProxyVecRNG();
   
   /** @brief Bare constructor - only for use with 'join' and 'split' methods */
   // JoiningProxyVecRng();

   Join( BaseVecRNGState_scalar* trackRngState[VecLen], int /*numVariatesExpected*/ );
   Join( SequentialRNG*          trackPrng[VecLen],     int /*numVariatesExpected*/ );
   Split();
   
private:
   // Location of track's RNG objects - in order to copy state back
   BaseRNGState_scalar* fSequentialTrackPrng[BackendT::VectorSize];
   // BaseRNGState_vector  fVectorRngState;
   BaseVecRNGType       fVectorRng;  //  Vector RNG which has 'sucked up' scalar state
};


// Implementation
template <typename BaseVecRNGType, typename BackendT>
  JoiningProxyVecRNG::
   Join( BaseRNGState* trackPrng[], int vecLen )
{
   assert( vecLen == BackendT::VectorSize );
   
   // Keep record of locations of per-track PRNG - in order to copy back the final state
   for( int i= 0; i < BackendT::VectorSize; ++i )
      fSequentialTrackPrng[i] = trackPrng[i];

   Initialize( fSequentialTrackPrng, BackendT::VectorSize );

   // NEED to copy ‘in’ the states of the RNGs 
   // Ideas for potential methods:
   // 1.) Copy the state of the array of 
   CopyState( fSequentialTrackPrng );

   //  Use an extension of the existing (VecRNG) methods
   //   SetState(State_t *state) { fState = state; }
   // e.g.
   //   SetState(State_s *state, int lane);

   // and the existing
   //  State_t* GetState() const { return fState; }

   for( int i= 0; i < BackendT::VectorSize; ++i )
      SetState( fSequentialTrackPrng[i].GetState(), i);
}

template <typename BaseVecRNGType, typename BackendT>
  JoiningProxyVecRNG::
   Split()
{
   for( int i= 0; i < BackendT::VectorSize; ++i )   
      fSequentialTrackPrng[i]
}

}
}
#endif
