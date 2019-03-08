
/**
 * Gamma Random ALgorithm from the implementation  from Marsaglia and
Tsang’s method (Marsaglia and Tsang, 2000).
The gamma random value greater than 1 can be extended
linearly to values of gamma less than one (Hung et al., 2015)
(second part of each algorithm)
 - this code is a part of VecRNG.h
 *
 */
// Vectorized algorithm
/**
 * Gamma: Gamma deviates - this code is a part of VecRNG.h
 */
// Gamma deviates
template <typename DerivedT>
template <typename BackendT>
VECCORE_ATT_HOST_DEVICE
typename BackendT::Double_v
VecRNG<DerivedT>::GammaScalar(typename BackendT::Double_v alpha,
                              typename BackendT::Double_v beta) {
  using Double_v = typename BackendT::Double_v;
  if (alpha > 1.0) {
    Double_v g = 0.0;
    Double_v d = alpha - 1.0 / 3.0;
    Double_v c = 1.0 / math::Sqrt(9.0 * d);
    Double_v x =
        static_cast<DerivedT *>(this)->template Gauss<BackendT>(0.0, 1.0);
    Double_v f = 1.0 + c * x;
    Double_v v = f * f * f;
    while (v <= 0.0) {
      x = static_cast<DerivedT *>(this)->template Gauss<BackendT>(0.0, 1.0);
      f = 1.0 + c * x;
      v = f * f * f;
    }
    Double_v u = static_cast<DerivedT *>(this)->template Uniform<BackendT>();
    bool flag;
    flag = false;
    do {
      if (u < (1.0 - 0.0331 * x * x * x * x)) {
        flag = true;
        g = d * v / beta;
        return g;

      } else {
        if (math::Log(u) < (0.5 * x * x + d * (1.0 - v + math::Log(v)))) {
          flag = true;
          g = d * v / beta;
          return g;
        } else {
          x = static_cast<DerivedT *>(this)->template Gauss<BackendT>(0.0, 1.0);
          f = 1.0 + c * x;
          v = f * f * f;
          while (v <= 0.0) {
            x = static_cast<DerivedT *>(this)->template Gauss<BackendT>(0.0,
                                                                        1.0);
            f = 1.0 + c * x;
            v = f * f * f;
          }
          u = static_cast<DerivedT *>(this)->template Uniform<BackendT>();
        }
      }
    } while (flag == false);
  } else {
    Double_v gamma = GammaScalar<BackendT>(alpha + 1, beta);
    Double_v u = static_cast<DerivedT *>(this)->template Uniform<BackendT>();
    Double_v e = 1.0 / alpha;
    gamma *= math::Pow(u, e);
    return gamma;
  }
}

template <typename DerivedT>
template <typename BackendT>
VECCORE_ATT_HOST_DEVICE
typename BackendT::Double_v
VecRNG<DerivedT>::Gamma(typename BackendT::Double_v alpha,
                        typename BackendT::Double_v beta) {
  using Double_v = typename BackendT::Double_v;
  int vsize = VectorSize<Double_v>();
  Double_v gamma(0.0);
  Mask<Double_v> gtOne(alpha > Double_v(1.0));
  if (MaskFull(gtOne)) { // all alpha > 1
    Double_v d = alpha - 1.0 / 3;
    Double_v c = 1.0 / math::Sqrt(9.0 * d);
    Double_v g =
        static_cast<DerivedT *>(this)->template Gauss<BackendT>(0.0, 1.0);
    Double_v z = c * g + 1.0;
    Double_v v = z * z * z;
    Mask<Double_v> good(v >= Double_v(0.0));

    if (MaskFull(
            good)) { // test for normal(0,1) < - sqrt(9 * alpha - 3) < -2.45
      Double_v u = static_cast<DerivedT *>(this)->template Uniform<BackendT>();
      Mask<Double_v> flag(u < (1 - 0.0331 * g * g * g * g));
      if (MaskFull(flag)) {
        // MaskedAssign(y, !done, 1.0 - x);
        gamma = d * v / beta;
        return gamma;
      } else {
        Double_v pe = math::Log(u);
        Double_v po = math::Log(v);
        Mask<Double_v> flag2(pe < (0.5 * g * g + d * (1.0 - v + po)));
        if (MaskFull(flag2)) {
          gamma = d * v / beta;
          return gamma;
        } else {
          for (int i = 0; i < vsize; ++i) {
            if (flag[i] == true) { // accepted
              gamma[i] = d[i] * v[i] / beta[i];
            } else {
              if (flag2[i] == true) { // accepted
                gamma[i] = d[i] * v[i] / beta[i];
              } else {
                //  gamma[i] = GammaScalar<ScalarBackend>(alpha[i],beta[i]);
                bool flag_s;
                flag_s = false;
                do {
                  // double f= static_cast<DerivedT *>(this)-> template
                  // Gauss<BackendT>(0.0,1.0); f=g[i];
                  g[i] = static_cast<DerivedT *>(this)
                             ->template Gauss<ScalarBackend>(0., 1.);
                  z[i] = c[i] * g[i] + 1.0;
                  v[i] = z[i] * z[i] * z[i];
                  if (v[i] >= 0.0) {
                    //  Mask<Double_v> good (v >= Double_v(0.0));
                    u[i] = static_cast<DerivedT *>(this)
                               ->template Uniform<ScalarBackend>();
                    //  pe[i]=math::Log(u[i]);
                    //  po[i]=math::Log(v[i]);
                    if (u[i] < (1.0 - 0.0331 * g[i] * g[i] * g[i] * g[i])) {
                      flag_s = true;
                    }
                  }
                } while (flag_s == false);
                gamma[i] = d[i] * v[i] / beta[i];
              }
            }
          }
          return gamma;
        }
      }
    } else {
      Double_v u = static_cast<DerivedT *>(this)->template Uniform<BackendT>();
      Double_v pe = math::Log(u);
      Double_v po = math::Log(v);
      for (int i = 0; i < vsize; ++i) {
        if (good[i] == true) { // accepted
          u[i] =
              static_cast<DerivedT *>(this)->template Uniform<ScalarBackend>();
          //
          if (u[i] < (1.0 - 0.0331 * g[i] * g[i] * g[i] * g[i])) {
            gamma[i] = d[i] * v[i] / beta[i];
          } else {
            if (pe[i] < (0.5 * g[i] * g[i] + d[i] * (1.0 - v[i] + po[i]))) {
              gamma[i] = d[i] * v[i] / beta[i];
            } else {
              bool flag_s;
              flag_s = false;
              do {
                g[i] = static_cast<DerivedT *>(this)
                           ->template Gauss<ScalarBackend>(0., 1.);
                z[i] = c[i] * g[i] + 1.0;
                v[i] = z[i] * z[i] * z[i];
                if (v[i] >= 0.0) {
                  //  Mask<Double_v> good (v >= Double_v(0.0));
                  u[i] = static_cast<DerivedT *>(this)
                             ->template Uniform<ScalarBackend>();
                  if (u[i] < (1.0 - 0.0331 * g[i] * g[i] * g[i] * g[i])) {
                    flag_s = true;
                  }
                }
              } while (flag_s == false);
              gamma[i] = d[i] * v[i] / beta[i];
            }
          }
        }
      }
      return gamma;
    }
  } else {                  // gama not all > 1
    if (MaskEmpty(gtOne)) { // 0 < alpha < 1
      Double_v gamma = Gamma<BackendT>(alpha + 1, beta);
      Double_v u = static_cast<DerivedT *>(this)->template Uniform<BackendT>();
      Double_v e = 1.0 / alpha;
      gamma *= math::Pow(u, e);
      return gamma;
    } else {
      // mixed cases for 0 < alpha < 1 and alpha > 1 - the scalar loop
      for (int i = 0; i < vsize; ++i) {
        gamma[i] = GammaScalar<ScalarBackend>(alpha[i], beta[i]);
      }
      // Gamma<ScalarBackend> may be specialized with the original
      // implementation
    }
  }
}
