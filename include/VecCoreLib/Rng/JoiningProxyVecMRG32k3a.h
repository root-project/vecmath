#ifndef JOINING_PROXY_VEC_MRG32k3A_H
#define JOINING_PROXY_VEC_MRG32k3A_H

#include "RngDefs.h"
#include "MRG32k3a.h"

namespace vecRng {
inline namespace VECRNG_IMPL_NAMESPACE {

template <typename BackendT>
   class JoiningProxyVecMRG32k3a :
     public MRG32k3a<BackendT>
     // public BaseVecRNGType /* <BackendT> */    // Type 1: derive & extend
     // public VecRNG< JoiningProxyVecMRG32k3a<BackendT> >  // Type 2: Use 'has'
{
public:
   using BaseVecRNGType=    MRG32k3a<VectorBackend>; 
   using BaseScalarRNGType= MRG32k3a<ScalarBackend>;
   
   // using State_t= typename RNG_traits<BaseVecRNGType> <BaseVecRNGType,BaseScalarRNGType,BackendT>> ::State_t;
   using BaseRNGState_scalar= typename RNG_traits<BaseVecRNGType> /*<ScalarBackend>>*/ ::State_t;
   using BaseRNGState_vector= typename RNG_traits<BaseScalarRNGType> /*<VectorBackend>>*/ ::State_t;
   using State_t= BaseRNGState_vector;

   /** @brief Constructor: 'create' empty instance - not associated with set of scalar states */
   JoiningProxyVecMRG32k3a();
   
   /** @brief Constructor: 'create' vector instance from a set of scalar states */
   JoiningProxyVecMRG32k3a( BaseScalarRNGType* trackPrng[], int numStates );
   // JoiningProxyVecMRG32k3a( BaseRNGState_scalar* trackPrng[], int vecLength );

   /** @brief Split apart this class - ensure that the state of its parts is fed back */
   //          - i.e. the state gets copied back to each 'track' scalar RNG state
   ~JoiningProxyVecMRG32k3a();
   
   /** @brief Bare constructor - only for use with 'join' and 'split' methods */
   // JoiningProxyVecRng();

   void Join( BaseScalarRNGType* trackRngState[], int numStates );
   // void Join( BaseScalarRNGState_scalar* trackRngState[], int numStates );   
   // Join( ScalarRNG*       trackPrng[VecLen],     int /*numVariatesExpected*/ );
   void Split();

   VECCORE_ATT_HOST_DEVICE   
   typename BackendT::Double_v Uniform();
   
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

/*   
  template <typename ReturnBackendType>
  VECCORE_ATT_HOST_DEVICE
  VECCORE_FORCE_INLINE
     typename ReturnBackendType::Double_v Kernel(State_t& state);
*/
   
private:
   // Location of track's RNG objects - in order to copy state back
   // BaseScalarRNGType*  fScalarTrackPrng[VectorSize<BackendT>];
   BaseScalarRNGType**  fScalarTrackPrng;
   bool                fFullState= false;
   // BaseRNGState_vector  fVectorRngState;
   // BaseVecRNGType       fVectorRng;  //  Vector RNG which has 'sucked up' scalar state
};

/*****
template <typename BaseVecRNGType,
          typename BaseScalarRNGType,
          typename BackendT>
1struct RNG_traits<JoiningProxyVecMRG32k3a<BaseVecRNGType, BaseScalarRNGType, BackendT>>
{
  struct StateJoined { BaseVecRNGType fVectorRng; }
  using State_t = StateJoined;
};
******/

// Implementation
template <typename BackendT>
void
  JoiningProxyVecMRG32k3a<BackendT>::
   Join( BaseScalarRNGType* trackPrng[], int numStates )
{
   assert( numStates == VectorSize<BackendT> );
        // If other parameters must be initialised ...

   for( int i= 0; i < BackendT::VectorSize; ++i ) {
      assert( trackPrng[i] != nullptr );

      SetState( fScalarTrackPrng[i].GetState(), i);
      // Keep record of locations of per-track PRNG - in order to copy back the final state         
      fScalarTrackPrng[i] = trackPrng[i];

   }
   fFullState= true;
}

//  ---------------------------================---------------------------
template <typename BackendT>
void JoiningProxyVecMRG32k3a<BackendT>::Split()
{
   for( int i= 0; i < BackendT::VectorSize; ++i ) {
      fScalarTrackPrng[i].SetState(
         MRG32k3a<BackendT>::GetState().fVectorRNG.GetPartOfState(i) ); // Sets full state
      fScalarTrackPrng[i] = nullptr;
   }
   fFullState= false;
}

//  ---------------------------================---------------------------
template <typename BackendT>
   JoiningProxyVecMRG32k3a<BackendT>::JoiningProxyVecMRG32k3a()
{
   // fScalarTrackPrng= new BaseScalarRNGType* [VectorSize<BackendT>] ;
   
   fFullState= false;  // Not
   MRG32k3a<BackendT>::Initialize(); //  fScalarTrackPrng, BackendT::VectorSize );
   for( int i= 0; i < BackendT::VectorSize; ++i )
      fScalarTrackPrng[i] = nullptr;   
}

//  ---------------------------================---------------------------
template <typename BackendT>
   JoiningProxyVecMRG32k3a<BackendT>::JoiningProxyVecMRG32k3a( BaseScalarRNGType* trackPrng[], int numStates )
   : JoiningProxyVecMRG32k3a()
{
   // Initialize();
   Join( trackPrng, numStates );
}

//  ---------------------------================---------------------------
template <typename BackendT>
   JoiningProxyVecMRG32k3a<BackendT>::~JoiningProxyVecMRG32k3a()
{
   if( fFullState ) Split();
   delete[] fScalarTrackPrng;
}

VECCORE_ATT_HOST_DEVICE
template <typename BackendT>
  typename BackendT::Double_v JoiningProxyVecMRG32k3a<BackendT>::Uniform()
{
   assert( fFullState );
   MRG32k3a<BackendT>::Uniform();
}

/***
// If redefinining the type completely
template <typename BackendT>
template <class ReturnTypeBackendT>   
typename ReturnTypeBackendT::Double_v
  JoiningProxyVecMRG32k3a<BackendT>::
   ::Kernel(State_t& state)
{
   BaseVecRNGType::Kernel( tstate.fBaseState );
}
****/

}
}
#endif
