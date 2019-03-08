#include "VecCore/Timer.h"

#include "VecMath/Rng/MRG32k3a.h"
#include "VecMath/Rng/Philox.h"
#include "VecMath/Rng/Threefry.h"

#ifdef RNGTEST_MKL
#include "mkl_vsl.h"
#include <omp.h>
#endif

#include "RngTest.h"

namespace vecrng {

// Scalar

double ScalarMRG32k3a(int nsample, double &result) {
  // Scalar MRG32k3a
  static vecRng::cxx::MRG32k3a<ScalarBackend> rng;
  rng.Initialize();

  static Timer<nanoseconds> timer;
  double elapsedTime = 0.;

  double sum = 0;

  timer.Start();

  for (int i = 0; i < nsample; ++i) {
    sum += rng.Uniform<ScalarBackend>();
  }

  elapsedTime = timer.Elapsed();
  result = sum;

  return elapsedTime;
}

double ScalarThreefry(int nsample, double &result) {
  // Scalar Threefry
  static vecRng::cxx::Threefry<ScalarBackend> rng;
  rng.Initialize();

  static Timer<nanoseconds> timer;
  double elapsedTime = 0.;

  double sum = 0;

  timer.Start();

  for (int i = 0; i < nsample; ++i) {
    sum += rng.Uniform<ScalarBackend>();
  }

  elapsedTime = timer.Elapsed();
  result = sum;

  return elapsedTime;
}

double ScalarPhilox(int nsample, double &result) {
  // Scalar Philox - use scalar method of VectorBackend (temporarily)

  static vecRng::cxx::Philox<ScalarBackend> rng;
  rng.Initialize();

  static Timer<nanoseconds> timer;
  double elapsedTime = 0.;

  double sum = 0;

  timer.Start();

  for (int i = 0; i < nsample; ++i) {
    sum += rng.Uniform<ScalarBackend>();
  }

  elapsedTime = timer.Elapsed();
  result = sum;

  vecRng::cxx::Philox<ScalarBackend> rng1 = rng;

  return elapsedTime;
}

// Vector

double VectorMRG32k3a(int nsample, double &result) {
  // Vector MRG32k3a
  using Double_v = typename VectorBackend::Double_v;
  int vsize = VectorSize<Double_v>();

  vecRng::cxx::MRG32k3a<VectorBackend> rng;
  rng.Initialize();

  static Timer<nanoseconds> timer;
  double elapsedTime = 0.;

  Double_v sum = 0.;

  timer.Start();

  for (int i = 0; i < nsample / vsize; ++i) {
    sum += rng.Uniform<VectorBackend>();
  }

  elapsedTime = timer.Elapsed();
  for (int i = 0; i < vsize; ++i)
    result += sum[i];

  return elapsedTime;
}

double VectorThreefry(int nsample, double &result) {
  // Vector Threefry
  using Double_v = typename VectorBackend::Double_v;
  int vsize = VectorSize<Double_v>();

  static vecRng::cxx::Threefry<VectorBackend> rng;
  rng.Initialize();

  static Timer<nanoseconds> timer;
  double elapsedTime = 0.;

  Double_v sum = 0;

  timer.Start();

  for (int i = 0; i < nsample / vsize; ++i) {
    sum += rng.Uniform<VectorBackend>();
  }

  elapsedTime = timer.Elapsed();
  for (int i = 0; i < vsize; ++i)
    result += sum[i];

  return elapsedTime;
}

double VectorPhilox(int nsample, double &result) {
  // Vector Philox
  using Double_v = typename VectorBackend::Double_v;
  int vsize = VectorSize<Double_v>();

  vecRng::cxx::Philox<VectorBackend> rng;
  rng.Initialize();

  static Timer<nanoseconds> timer;
  double elapsedTime = 0.;

  Double_v sum = 0;

  timer.Start();

  for (int i = 0; i < nsample / vsize; ++i) {
    sum += rng.Uniform<VectorBackend>();
  }

  elapsedTime = timer.Elapsed();
  for (int i = 0; i < vsize; ++i)
    result += sum[i];

  return elapsedTime;
}

double StateMRG32k3a(int nsample, double &result) {
  // Scalar MRG32k3a
  static vecRng::cxx::MRG32k3a<ScalarBackend> rng;

  int theNBlocks = 26;
  int theNThreads = 192;

  vecRng::MRG32k3a<ScalarBackend>::State_t *hstates =
      (vecRng::MRG32k3a<ScalarBackend>::State_t *)malloc(
          theNBlocks * theNThreads *
          sizeof(vecRng::MRG32k3a<ScalarBackend>::State_t));
  rng.Initialize(hstates, theNBlocks * theNThreads);

  static Timer<nanoseconds> timer;
  double elapsedTime = 0.;

  double sum = 0;

  timer.Start();

  for (int i = 0; i < theNBlocks * theNThreads; ++i) {
    for (int j = 0; j < nsample / (theNBlocks * theNThreads); ++j) {
      sum += rng.Uniform<ScalarBackend>(&hstates[i]);
    }
  }

  elapsedTime = timer.Elapsed();
  result = sum;

  free(hstates);

  return elapsedTime;
}

double StateThreefry(int nsample, double &result) {
  // Scalar Threefry
  static vecRng::cxx::Threefry<ScalarBackend> rng;

  int theNBlocks = 26;
  int theNThreads = 192;

  vecRng::Threefry<ScalarBackend>::State_t *hstates =
      (vecRng::Threefry<ScalarBackend>::State_t *)malloc(
          theNBlocks * theNThreads *
          sizeof(vecRng::Threefry<ScalarBackend>::State_t));
  rng.Initialize(hstates, theNBlocks * theNThreads);

  static Timer<nanoseconds> timer;
  double elapsedTime = 0.;

  double sum = 0;

  timer.Start();

  for (int i = 0; i < theNBlocks * theNThreads; ++i) {
    for (int j = 0; j < nsample / (theNBlocks * theNThreads); ++j) {
      sum += rng.Uniform<ScalarBackend>(&hstates[i]);
    }
  }

  elapsedTime = timer.Elapsed();
  result = sum;

  free(hstates);

  return elapsedTime;
}

double StatePhilox(int nsample, double &result) {
  // Scalar Philox
  //  static vecrng::cxx::VecRNG<Philox<ScalarBackend>, ScalarBackend,
  //  Philox_t<ScalarBackend>> rng;
  static vecRng::cxx::Philox<ScalarBackend> rng;

  int theNBlocks = 26;
  int theNThreads = 192;

  vecRng::Philox<ScalarBackend>::State_t *hstates =
      (vecRng::Philox<ScalarBackend>::State_t *)malloc(
          theNBlocks * theNThreads *
          sizeof(vecRng::Philox<ScalarBackend>::State_t));
  rng.Initialize(hstates, theNBlocks * theNThreads);

  static Timer<nanoseconds> timer;
  double elapsedTime = 0.;

  double sum = 0;

  timer.Start();

  for (int i = 0; i < theNBlocks * theNThreads; ++i) {
    for (int j = 0; j < nsample / (theNBlocks * theNThreads); ++j) {
      sum += rng.Uniform<ScalarBackend>(&hstates[i]);
    }
  }

  elapsedTime = timer.Elapsed();
  result = sum;

  free(hstates);

  return elapsedTime;
}

// Intel VMK/VSL
#ifdef RNGTEST_MKL

double VSLRngTest(int nsample, double &result, int rngType) {
  static Timer<nanoseconds> timer;
  double elapsedTime = 0.;

  // buffer size
  const int N = 32;

  double r[N];
  double a = 0.0;
  double b = 1.0;

  double sum = 0;

  // Initialize
  unsigned int seed = 7777777;
  VSLStreamStatePtr stream;
  vslNewStream(&stream, rngType, (MKL_INT)seed);

  timer.Start();

  // Call RNG
  for (int i = 0; i < nsample / N; ++i) {
    vdRngUniform(VSL_RNG_METHOD_UNIFORM_STD_ACCURATE, stream, N, r, a, b);
    for (int j = 0; j < N; ++j)
      sum += r[j];
  }
  elapsedTime = timer.Elapsed();

  result = sum;

  // Deinitialize
  vslDeleteStream(&stream);

  return elapsedTime;
}

#endif

} // end namespace vecrng
