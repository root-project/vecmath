#ifndef JOINING_PROXY_VEC_RNG_H
#define JOINING_PROXY_VEC_RNG_H

#include "RngDefs.h"
#include "VecRNG.h"

namespace vecRng {
inline namespace VECRNG_IMPL_NAMESPACE {

template <typename BaseVecRNGType,
          typename BaseScalarRNGType,
          typename BackendT>
   class JoiningProxyVecRNG :
     // public BaseVecRNGType /* <BackendT> */    // Type 1: derive & extend
     public VecRNG< JoiningProxyVecRNG<BaseVecRNGType, BaseScalarRNGType, BackendT> >  // Type 2: Use 'has'
{
public:
   // using State_t= typename RNG_traits<BaseVecRNGType> <BaseVecRNGType,BaseScalarRNGType,BackendT>> ::State_t;
   using BaseRNGState_scalar= typename RNG_traits<BaseVecRNGType> /*<ScalarBackend>>*/ ::State_t;
   using BaseRNGState_vector= typename RNG_traits<BaseScalarRNGType> /*<VectorBackend>>*/ ::State_t;
   using State_t= BaseRNGState_vector;
   
   /** @brief Constructor: 'create' vector instance from a set of sequential states */
   JoiningProxyVecRNG( BaseScalarRNGType* trackPrng[], int numStates );
   // JoiningProxyVecRNG( BaseRNGState_scalar* trackPrng[], int vecLength );

   /** @brief Split apart this class - ensure that the state of its parts is fed back */
   //          - i.e. the state gets copied back to each 'track' sequential RNG state
   ~JoiningProxyVecRNG();
   
   /** @brief Bare constructor - only for use with 'join' and 'split' methods */
   // JoiningProxyVecRng();

   void Join( BaseScalarRNGType* trackRngState[], int numStates );
   // void Join( BaseScalarRNGState_scalar* trackRngState[], int numStates );   
   // Join( SequentialRNG*       trackPrng[VecLen],     int /*numVariatesExpected*/ );
   void Split();

/*****   
  // Initialization for SIMD
  VECCORE_ATT_HOST
  void Initialize()
  { static_cast<DerivedT *>(this)->template Initialize(); }

  // Initialization with a unique stream number
  VECCORE_ATT_HOST
  void Initialize(long streamId)
  { static_cast<DerivedT *>(this)->template Initialize(streamId); }

  // Initialization for SIMT
  VECCORE_ATT_HOST
  void Initialize(State_t *states, unsigned int nthreads)
  { static_cast<DerivedT *>(this)->template Initialize(states, nthreads); }
 *****/
   
  // ------------------------------------------------------------------------------------    
  // Expected methods for VecRNG ... needed in 'has' design
  // 

  template <typename ReturnBackendType>
  VECCORE_ATT_HOST_DEVICE
  VECCORE_FORCE_INLINE
     typename ReturnBackendType::Double_v Kernel(State_t& state);

/*********************  
  // ------------------------------------------------------------------------------------  
  // Methods provided by VecRNG ... no need in 'has' design 
   
  // Returns an array of random numbers of the type BackendT::Double_v
  template <typename BackendT>
  VECCORE_ATT_HOST_DEVICE
  void Array(const size_t nsize, typename BackendT::Double_v *array)
  { fVectorRng.Array<BackendT>; }

  // Flat distribution in [min,max)
  template <typename BackendT>
  VECCORE_ATT_HOST_DEVICE
  typename BackendT::Double_v Flat(typename BackendT::Double_v min,
                                   typename BackendT::Double_v max)
  { fVectorRng.Flat<BackendT>(min,max); }

  // Flat distribution in [min,max] with a state
  template <typename BackendT>
  VECCORE_ATT_HOST_DEVICE
  typename BackendT::Double_v Flat(State_t *state,
                                   typename BackendT::Double_v min,
                                   typename BackendT::Double_v max)  
  { return min+(max-min)*static_cast<DerivedT *>(this)-> template Uniform<BackendT>(state); }
 *******************/

private:
   // Location of track's RNG objects - in order to copy state back
   // BaseRNGState_scalar* fSequentialTrackPrng[BackendT::VectorSize];
   BaseScalarRNGType*  fSequentialTrackPrng[BackendT::VectorSize];
   
   // BaseRNGState_vector  fVectorRngState;
   // BaseVecRNGType       fVectorRng;  //  Vector RNG which has 'sucked up' scalar state
};

template <typename BaseVecRNGType,
          typename BaseScalarRNGType,
          typename BackendT>
struct RNG_traits<JoiningProxyVecRNG<BaseVecRNGType, BaseScalarRNGType, BackendT>>
{
  // fCg[row=MRG::vsize][column=VectorSize<Double_v>()]
  struct StateJoined { BaseVecRNGType fVectorRng; }
  using State_t = StateJoined;
};

// Implementation
template <typename BaseVecRNGType,
          typename BaseScalarRNGType,   
          typename BackendT>
void
  JoiningProxyVecRNG<BaseVecRNGType, BaseScalarRNGType, BackendT>::
   Join( BaseScalarRNGType* trackPrng[], int numStates )
{
   assert( numStates == VectorSize<BackendT> ); //  BackendT::VectorSize );
   
   // Keep record of locations of per-track PRNG - in order to copy back the final state
   for( int i= 0; i < BackendT::VectorSize; ++i )
      fSequentialTrackPrng[i] = trackPrng[i];

   VecRng<BackendT>::Initialize(); // fSequentialTrackPrng, BackendT::VectorSize );

   // NEED to copy ‘in’ the states of the RNGs 
   // Ideas for potential methods:
   // 1.) Copy the state of the array of 
   // CopyState( fSequentialTrackPrng );

   //  Use an extension of the existing (VecRNG) methods
   //   SetState(State_t *state) { fState = state; }
   // e.g.
   //   SetState(State_s *state, int lane);

   // and the existing
   //  State_t* GetState() const { return fState; }

   for( int i= 0; i < BackendT::VectorSize; ++i )
      SetState( fSequentialTrackPrng[i].GetState(), i);
}

//  ---------------------------================---------------------------

template <typename BaseVecRNGType,
          typename BaseScalarRNGType,   
          typename BackendT>
void JoiningProxyVecRNG<BaseVecRNGType, BaseScalarRNGType, BackendT>::Split()
{
   for( int i= 0; i < BackendT::VectorSize; ++i )   
      fSequentialTrackPrng[i].SetState( VecRNG<BackendT>::GetState().fVectorRNG.GetPartOfState(i) ); // Sets full state
                                        
}

//  ---------------------------================---------------------------

template <typename BaseVecRNGType,
          typename BaseScalarRNGType,   
          typename BackendT>
template <class ReturnTypeBackendT>   
typename ReturnTypeBackendT::Double_v
  JoiningProxyVecRNG<BaseVecRNGType, BaseScalarRNGType, BackendT>::
   ::Kernel(State_t& state)
{
   BaseVecRNGType::Kernel( tstate.fBaseState );
}

}
}
#endif
