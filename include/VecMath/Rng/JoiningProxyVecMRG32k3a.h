#ifndef JOINING_PROXY_VEC_MRG32k3A_H
#define JOINING_PROXY_VEC_MRG32k3A_H

#include "RngDefs.h"
#include "MRG32k3a.h"

namespace vecRng {
inline namespace VECRNG_IMPL_NAMESPACE {

template <typename BackendT>
   class JoiningProxyVecMRG32k3a :
        public MRG32k3a<BackendT>
{
public:
   using BaseVecRNGType=    MRG32k3a<VectorBackend>; 
   using BaseScalarRNGType= MRG32k3a<ScalarBackend>;
   
   using BaseRNGState_scalar= typename RNG_traits<BaseScalarRNGType> /*<ScalarBackend>>*/ ::State_t;
   using BaseRNGState_vector= typename RNG_traits<BaseVecRNGType> /*<VectorBackend>>*/ ::State_t;
   using State_t= BaseRNGState_vector;

   /** @brief Constructor: 'create' empty instance - not associated with set of scalar states */
   JoiningProxyVecMRG32k3a();
   
   /** @brief Constructor: 'create' vector instance from a set of scalar states */
   JoiningProxyVecMRG32k3a( BaseScalarRNGType* trackPrng[], int numStates );

   /** @brief Split apart this class - ensure that the state of its parts is fed back */
   //          - i.e. the state gets copied back to each 'track' scalar RNG state
   ~JoiningProxyVecMRG32k3a();

   /** @brief Create 'vector' state by copying the array of per-track states */
   void Join( BaseScalarRNGType* trackRngState[], int numStates );

   /** @brief De-couple 'vector' state - copy back the per-track states into original arrays */
   void Split();

   VECCORE_ATT_HOST_DEVICE   
   typename BackendT::Double_v Uniform();

   void PrintComponents() const;
   
protected:
   static typename BackendT::Double_v  dbVec;
   static constexpr int vecSize= VectorSize<decltype(dbVec)>();
   
private:
   // Location of track's RNG objects - in order to copy state back
   // BaseScalarRNGType*  fScalarTrackPrng[VectorSize<BackendT>];
   BaseScalarRNGType**  fScalarTrackPrng;
   bool                fFullState= false;
};

// Implementation
template <typename BackendT>
void
  JoiningProxyVecMRG32k3a<BackendT>::
   Join( BaseScalarRNGType* trackPrng[], int numStates )
{
   //  Recommendation to use VectorSize<decltype(x)>(); with x a variable of relevant type;
   typename BackendT::Double_v  dbVec;
   constexpr int vecSize= VectorSize<decltype(dbVec)>(); // = VectorSize<BackendT>();
   assert( numStates == vecSize );
   assert( fScalarTrackPrng );
        // If other parameters must be initialised ...

   // std::cout << "Locations of per-track PRNG: " << std::endl;
   
   for( int i= 0; i < vecSize; ++i ) {
      assert( trackPrng[i] != nullptr );
      
      MRG32k3a<BackendT>::SetPartOfState( *trackPrng[i]->GetState(), i);
      // Keep record of locations of per-track PRNG - in order to copy back the final state
      fScalarTrackPrng[i] = trackPrng[i];

      // std::cout << " [ " << i << " ] : " << trackPrng[i] << " stored: " << fScalarTrackPrng[i] << std::endl;
   }
   fFullState= true;
}

//  ---------------------------================---------------------------
template <typename BackendT>
void JoiningProxyVecMRG32k3a<BackendT>::Split()
{
   // std::cout << "Split called for object " << this << " State addr = " << this->fState
   //          << " vecSize= " << vecSize << std::endl;
   // std::cout << "Split called - vecSize = " << vecSize << std::endl;
   // PrintComponents(); 
   for( int i= 0; i < vecSize; ++i ) {
      // std::cout << " [ " << i << " ] : " << fScalarTrackPrng[i] << std::endl;
      
      auto partState = MRG32k3a<BackendT>::GetPartOfState(i); // Sets full state
      // Copy the slice's state back into the scalar 
      // *fScalarTrackPrng[i] = MRG32k3a<ScalarBackend>(partState);
      fScalarTrackPrng[i]->CopyState(partState);
      fScalarTrackPrng[i] = nullptr;
   }
   fFullState= false;
}

//  ---------------------------================---------------------------

template <typename BackendT>
void JoiningProxyVecMRG32k3a<BackendT>::PrintComponents() const
{
   const char* typeNameStr = "JoiningProxyVecMRG32k3a< >";
   std::cout << std::endl;
   std::cout << "== Printing of JoiningProxVecMRG32k3a PrintComponents():  " << std::endl;
   std::cout << "===================================================================" << std::endl;
   std::cout << "Dump of " << typeNameStr << " < Backend with vecSize= " << vecSize
             << " >  called for object " << this << " State addr = " << this->fState
             << std::endl;
   std::cout << " fullState = " << fFullState << std::endl;
   std::cout << " Contents of scalar RNG addresses: " << std::endl;
   for( int i= 0; i < vecSize; ++i ) {
      std::cout << " Addr [ " << i << " ] : " << fScalarTrackPrng[i] << std::endl;
   }
   std::cout << "===================================================================" << std::endl;   
}

//  ---------------------------================---------------------------
template <typename BackendT>
   JoiningProxyVecMRG32k3a<BackendT>::JoiningProxyVecMRG32k3a()
{
   fScalarTrackPrng= new BaseScalarRNGType* [vecSize];
   
   fFullState= false;  // Not
   constexpr int vecSize= VectorSize<BackendT>();
   MRG32k3a<BackendT>::Initialize(); //  fScalarTrackPrng, BackendT::VectorSize );
   for( int i= 0; i < vecSize; ++i )
      fScalarTrackPrng[i] = nullptr;   
}

//  ---------------------------================---------------------------
template <typename BackendT>
   JoiningProxyVecMRG32k3a<BackendT>::JoiningProxyVecMRG32k3a( BaseScalarRNGType* trackPrng[], int numStates )
   : JoiningProxyVecMRG32k3a<BackendT>()
{
   MRG32k3a<BackendT>::Initialize();
   fFullState= false;
   
   fScalarTrackPrng= new BaseScalarRNGType* [vecSize];   
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
   return this->MRG32k3a<BackendT>::template Uniform<BackendT>();
}

}
}
#endif
