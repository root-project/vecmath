/**
 * Gauss: Gaussian deviates - this code is a part of VecRNG.h
 *
 */

// Gaussian deviates with a state
template <typename DerivedT>
template <typename ReturnTypeBackendT>
VECCORE_ATT_HOST_DEVICE
typename ReturnTypeBackendT::Double_v
VecRNG<DerivedT>::Gauss(State_t *state,
                        typename ReturnTypeBackendT::Double_v mean,
                        typename ReturnTypeBackendT::Double_v sigma) {
  // Gauss with a state
  using Double_v = typename ReturnTypeBackendT::Double_v;

  Double_v u1 =
      static_cast<DerivedT *>(this)->template Uniform<ReturnTypeBackendT>(
          state);
  Double_v u2 =
      static_cast<DerivedT *>(this)->template Uniform<ReturnTypeBackendT>(
          state) *
      (2.0 * M_PI);
  Double_v normal = math::Sqrt(-2.0 * math::Log(u1)) * math::Cos(u2);
  return mean + sigma * normal;
}

// Gaussian deviates with (mean, sigma):
// 1/(2*pi*sigma^2)*exp[-(x-mean)^2/sigma^2]
template <typename DerivedT>
template <typename ReturnTypeBackendT>
VECCORE_ATT_HOST_DEVICE
typename ReturnTypeBackendT::Double_v
VecRNG<DerivedT>::Gauss(typename ReturnTypeBackendT::Double_v mean,
                        typename ReturnTypeBackendT::Double_v sigma) {
  // Using Box/Muller - use just one
  // normal1 = sqrt(-2*log(u1))*cos(2*pi*u2)
  // normal2 = sqrt(-2*log(u1))*sin(2*pi*u2)

  using Double_v = typename ReturnTypeBackendT::Double_v;

  Double_v u1 =
      static_cast<DerivedT *>(this)->template Uniform<ReturnTypeBackendT>();
  Double_v u2 =
      static_cast<DerivedT *>(this)->template Uniform<ReturnTypeBackendT>() *
      (2.0 * M_PI);
  Double_v normal = math::Sqrt(-2.0 * math::Log(u1)) * math::Cos(u2);
  return mean + sigma * normal;
}
