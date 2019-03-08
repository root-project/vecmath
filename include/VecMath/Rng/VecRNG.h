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

template <typename DerivedT> class RNG_traits;

template <typename DerivedT> class VecRNG {

protected:
  // Use *this to access data members in the derived class
  using State_t = typename RNG_traits<DerivedT>::State_t;
  State_t *fState;

public:
  VECCORE_ATT_HOST_DEVICE
#ifndef VECCORE_CUDA
  VecRNG() {
    fState = (State_t *)AlignedAlloc(VECCORE_SIMD_ALIGN, sizeof(State_t));
  }
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
  void Initialize() { static_cast<DerivedT *>(this)->template Initialize(); }

  // Initialization with a unique stream number
  VECCORE_ATT_HOST
  void Initialize(long streamId) {
    static_cast<DerivedT *>(this)->template Initialize(streamId);
  }

  // Initialization for SIMT
  VECCORE_ATT_HOST
  void Initialize(State_t *states, unsigned int nthreads) {
    static_cast<DerivedT *>(this)->template Initialize(states, nthreads);
  }

  // Return BackendT::Double_v of random numbers in [0,1)
  template <typename BackendT> VECCORE_ATT_HOST_DEVICE typename BackendT::Double_v Uniform() {
    return static_cast<DerivedT *>(this)->template Kernel<BackendT>(
        *this->fState);
  }

  // Generate random numbers based on a given state
  template <typename BackendT>
  VECCORE_ATT_HOST_DEVICE
  typename BackendT::Double_v Uniform(State_t *state) {
    return static_cast<DerivedT *>(this)->template Kernel<BackendT>(*state);
  }

  VECCORE_ATT_HOST
  void PrintState() const { static_cast<DerivedT *>(this)->PrintState(); }

  // Auxiliary methods

  VECCORE_ATT_HOST_DEVICE
  void SetState(State_t *state) { fState = state; }

  VECCORE_ATT_HOST_DEVICE
  State_t *GetState() const { return fState; }

  VECCORE_ATT_HOST_DEVICE
  State_t const &GetStateRef() const { return *fState; }

  // Common methods

  // Return Index_v<ypename BackendT::Double_v> of random index numbers in
  // (min,max]
  template <typename BackendT>
  VECCORE_ATT_HOST_DEVICE
  Index_v<typename BackendT::Double_v>
  UniformIndex(Index_v<typename BackendT::Double_v> min = 0,
               Index_v<typename BackendT::Double_v> max = INT32_MAX) {
    return min +
           (max - min) *
               static_cast<DerivedT *>(this)->template Uniform<BackendT>();
  }

  // UniformIndex - specialization for scalar
  VECCORE_ATT_HOST_DEVICE
  Index_v<double> UniformIndex(Index_v<double> min = 0,
                               Index_v<double> max = UINT64_MAX) {
    return min +
           (max - min) *
               static_cast<DerivedT *>(this)->template Uniform<ScalarBackend>();
  }

  // Return Index_v<ypename BackendT::Double_v> of random numbers in (min,max]
  // with a state
  template <typename BackendT>
  VECCORE_ATT_HOST_DEVICE
  Index_v<typename BackendT::Double_v>
  UniformIndex(State_t *state, Index_v<typename BackendT::Double_v> min = 0,
               Index_v<typename BackendT::Double_v> max = INT32_MAX) {
    return min +
           (max - min) *
               static_cast<DerivedT *>(this)->template Uniform<BackendT>(state);
  }

  // UniformIndex with a status - specialization for scalar
  VECCORE_ATT_HOST_DEVICE
  Index_v<double> UniformIndex(State_t *state, Index_v<double> min = 0,
                               Index_v<double> max = UINT64_MAX) {
    return min +
           (max - min) *
               static_cast<DerivedT *>(this)->template Uniform<ScalarBackend>(
                   state);
  }

  // Returns an array of random numbers of the type BackendT::Double_v
  template <typename BackendT>
  VECCORE_ATT_HOST_DEVICE
  void Array(const size_t nsize, typename BackendT::Double_v *array);

  // Flat distribution in [min,max)
  template <typename BackendT>
  VECCORE_ATT_HOST_DEVICE
  typename BackendT::Double_v Flat(typename BackendT::Double_v min,
                                   typename BackendT::Double_v max) {
    return min +
           (max - min) *
               static_cast<DerivedT *>(this)->template Uniform<BackendT>();
  }

  // Flat distribution in [min,max] with a state
  template <typename BackendT>
  VECCORE_ATT_HOST_DEVICE
  typename BackendT::Double_v Flat(State_t *state,
                                   typename BackendT::Double_v min,
                                   typename BackendT::Double_v max) {
    return min +
           (max - min) *
               static_cast<DerivedT *>(this)->template Uniform<BackendT>(state);
  }

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
  // Gamma
  template <typename BackendT>
  VECCORE_ATT_HOST_DEVICE
  typename BackendT::Double_v Gamma(typename BackendT::Double_v alpha,
                                    typename BackendT::Double_v beta);

  // Gamma Scalar
  template <typename BackendT>
  VECCORE_ATT_HOST_DEVICE
  typename BackendT::Double_v GammaScalar(typename BackendT::Double_v alpha,
                                          typename BackendT::Double_v beta);

  // add Gamma with a state
};

// Implementation
// Common Methods
// Returns an array of random numbers of BackendT::Double_v
template <typename DerivedT>
template <typename BackendT>
VECCORE_ATT_HOST_DEVICE
void VecRNG<DerivedT>::Array(const size_t nsize,
                             typename BackendT::Double_v *array) {
  using Double_v = typename BackendT::Double_v;
  for (size_t i = 0; i < nsize; ++i) {
    Double_v u01 = static_cast<DerivedT *>(this)->template Uniform<BackendT>();
    array[i] = u01;
  }
}

// Exponential deviates: exp(-x/tau)
template <typename DerivedT>
template <typename BackendT>
VECCORE_ATT_HOST_DEVICE
typename BackendT::Double_v
VecRNG<DerivedT>::Exp(typename BackendT::Double_v tau) {
  using Double_v = typename BackendT::Double_v;

  Double_v u01 = static_cast<DerivedT *>(this)->template Uniform<BackendT>();
  //@syj: check for zero
  return -tau * math::Log(u01);
}

// Exponential deviates with a state
template <typename DerivedT>
template <typename BackendT>
VECCORE_ATT_HOST_DEVICE
typename BackendT::Double_v
VecRNG<DerivedT>::Exp(State_t *state, typename BackendT::Double_v tau) {
  // Exp with a state
  using Double_v = typename BackendT::Double_v;
  Double_v u01 =
      static_cast<DerivedT *>(this)->template Uniform<BackendT>(state);
  return -tau * math::Log(u01);
}
// INclude Gamma and Gauss PDF's generators
#include "Gamma.h"
#include "Gauss.h"

} // namespace VECRNG_IMPL_NAMESPACE
} // end namespace vecRng

#endif
