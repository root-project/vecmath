#ifndef VECRNG_VECRNG_H
#define VECRNG_VECRNG_H 1

/**
 * VecRNG: The base class of SIMD/SIMT random number generators
 *
 * Requirements :
 * 1) DerivedT  : A pseudo-random number generator with multiple streams
 * 2) State_t   : A templated struct of DerivedT states
 */

#include "RngDefs.h"

namespace vecRng {
inline namespace VECRNG_IMPL_NAMESPACE {

template <typename DerivedT>
struct RNG_traits;

/****
template <typename Derived>
struct RNG_traits<Derived> {
};
****/

template <typename DerivedT>
class VecRNG {

protected:
  // Use *this to access data members in the derived class
  using State_t = typename RNG_traits<DerivedT>::State_t;
  State_t *fState;

public:

  VECCORE_ATT_HOST_DEVICE
#ifndef VECCORE_CUDA
  VecRNG() { fState = (State_t *)AlignedAlloc(VECCORE_SIMD_ALIGN, sizeof(State_t)); }
#else
  VecRNG() { fState = new State_t; }
#endif

  // Dummy Constructor for SIMT
  VECCORE_ATT_HOST_DEVICE
  VecRNG(State_t *devState) { fState = devState; }

  VECCORE_ATT_HOST_DEVICE
#ifndef VECCORE_CUDA
  ~VecRNG() { AlignedFree(fState); }
#else
  ~VecRNG() { delete fState; }
#endif

  VECCORE_ATT_HOST_DEVICE
  VecRNG(const VecRNG &rng) = default;

  // Static interfaces (Required methods)

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

  // Return BackendT::Double_v of random numbers in [0,1)
  template <typename BackendT>
  VECCORE_ATT_HOST_DEVICE
  typename BackendT::Double_v Uniform()
  { return static_cast<DerivedT *>(this)->template Kernel<BackendT>(*this->fState); }

  // Generate random numbers based on a given state
  template <typename BackendT>
  VECCORE_ATT_HOST_DEVICE
  typename BackendT::Double_v Uniform(State_t *state)
  { return static_cast<DerivedT *>(this)->template Kernel<BackendT>(*state); }

  VECCORE_ATT_HOST
  void PrintState() const { static_cast<DerivedT *>(this)->PrintState(); }

  VECCORE_ATT_HOST_DEVICE
  void SetState(State_t *state) { fState = state; }

  VECCORE_ATT_HOST_DEVICE
  State_t* GetState() const { return fState; }

// #define DEFINE_PER_ELEMENT_GET_SET  1
#ifdef DEFINE_PER_ELEMENT_GET_SET  
  // Type for scalar VecRNG state
  // using ScalarBackend = vecCore::backend::Scalar;
  // using ScalarRngType = typename VecRNG<ScalarBackend>;
  // using State_scalar = typename RNG_traits<ScalarRngType>::State_t;    
  // using State_scalar = typename RNG_traits<DerivedT<ScalarBackend> >::State_t;
  // using State_v = typename RNG_traits<DerivedT<VectorBackend> >::State_t;

  // From VecCore backend/Interface.h: 
  // template <typename T>
  //  using Scalar = typename TypeTraits<T>::ScalarType;
  
  VECCORE_ATT_HOST_DEVICE
  void SetPartOfState( typename RNG_traits<VecRNG<ScalarBackend>>::State_t *state, // State_scalar
                       int i)
  { static_cast<DerivedT *>(this)->template SetState<BackendT>(state, i); }
     
  VECCORE_ATT_HOST_DEVICE
  typename RNG_traits<VecRNG<ScalarBackend>>::State_t * // State_scalar*
     GetPartOfState(int i) const
  { return static_cast<DerivedT *>(this)->template GetState<BackendT>(i); }
#endif
  // Auxiliary methods
  
  //Common methods

  // Returns an array of random numbers of the type BackendT::Double_v
  template <typename BackendT>
  VECCORE_ATT_HOST_DEVICE
  void Array(const size_t nsize, typename BackendT::Double_v *array);

  // Flat distribution in [min,max)
  template <typename BackendT>
  VECCORE_ATT_HOST_DEVICE
  typename BackendT::Double_v Flat(typename BackendT::Double_v min,
                                   typename BackendT::Double_v max)
  { return min+(max-min)*static_cast<DerivedT *>(this)-> template Uniform<BackendT>(); }

  // Flat distribution in [min,max] with a state
  template <typename BackendT>
  VECCORE_ATT_HOST_DEVICE
  typename BackendT::Double_v Flat(State_t *state,
                                   typename BackendT::Double_v min,
                                   typename BackendT::Double_v max)
  { return min+(max-min)*static_cast<DerivedT *>(this)-> template Uniform<BackendT>(state); }

  // Exponential deviates: exp(-x/tau)
  template <typename BackendT>
  VECCORE_ATT_HOST_DEVICE
  typename BackendT::Double_v Exp(typename BackendT::Double_v tau);

  // Exponential deviates with a state
  template <typename BackendT>
  VECCORE_ATT_HOST_DEVICE
  typename BackendT::Double_v Exp(State_t *state,
                                  typename BackendT::Double_v tau);

  // Gaussin deviates: 1/(2*pi*sigma^2)*exp[-(x-mean)^2/sigma^2]
  template <typename BackendT>
  VECCORE_ATT_HOST_DEVICE
  typename BackendT::Double_v Gauss(typename BackendT::Double_v mean,
                                    typename BackendT::Double_v sigma);

  // Gaussin deviates with a state
  template <typename BackendT>
  VECCORE_ATT_HOST_DEVICE
  typename BackendT::Double_v Gauss(State_t *state,
                                    typename BackendT::Double_v mean,
                                    typename BackendT::Double_v sigma);

  // @syj: add methods to generate other random distributions
  // (Binomial, Chi-Square, Gamma, Poisson, Landau and etc)

};

// Implementation

// Common Methods

// Returns an array of random numbers of BackendT::Double_v
template <typename DerivedT>
template <typename BackendT>
VECCORE_ATT_HOST_DEVICE void
VecRNG<DerivedT>::Array(const size_t nsize, typename BackendT::Double_v *array)
{
  using Double_v = typename BackendT::Double_v;
  for (size_t i = 0; i < nsize ; ++i) {
    Double_v u01 = static_cast<DerivedT *>(this)-> template Uniform<BackendT>();
    array[i] = u01;
  }
}

// Exponential deviates: exp(-x/tau)
template <typename DerivedT>
template <typename BackendT>
VECCORE_ATT_HOST_DEVICE typename BackendT::Double_v
VecRNG<DerivedT>::Exp(typename BackendT::Double_v tau)
{
  using Double_v = typename BackendT::Double_v;

  Double_v u01 = static_cast<DerivedT *>(this)-> template Uniform<BackendT>();
  //@syj: check for zero
  return -tau*math::Log(u01);
}

// Exponential deviates with a state
template <typename DerivedT>
template <typename BackendT>
VECCORE_ATT_HOST_DEVICE typename BackendT::Double_v
VecRNG<DerivedT>::Exp(State_t *state, typename BackendT::Double_v tau)
{
  // Exp with a state
  using Double_v = typename BackendT::Double_v;

  Double_v u01 = static_cast<DerivedT *>(this)-> template Uniform<BackendT>(state);
  return -tau*math::Log(u01);
}

// Gaussian deviates
#include "Gauss.h"

} // end namespace impl
} // end namespace vecRng

#endif
