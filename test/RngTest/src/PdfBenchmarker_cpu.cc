#include "VecCore/Timer.h"

#include "VecCoreLib/Rng/MRG32k3a.h"

#include "RngTest.h"

namespace vecrng {

// Scalar

double ScalarMRG32k3aExp(int nsample, double& result)
{
  // Scalar MRG32k3a Exp
  static vecRng::cxx::MRG32k3a<ScalarBackend> rng;
  rng.Initialize();

  static Timer<nanoseconds> timer;
  double elapsedTime = 0.;

  double sum = 0;

  timer.Start();

  for (int i = 0; i < nsample ; ++i) {
    sum += rng.Exp<ScalarBackend>(1.0);
  }

  elapsedTime = timer.Elapsed();
  result = sum;

  return elapsedTime;
}

double ScalarMRG32k3aNormal(int nsample, double& result)
{
  // Scalar MRG32k3a Normal
  static vecRng::cxx::MRG32k3a<ScalarBackend> rng;
  rng.Initialize();

  static Timer<nanoseconds> timer;
  double elapsedTime = 0.;

  double sum = 0;

  timer.Start();

  for (int i = 0; i < nsample ; ++i) {
    sum += rng.Gauss<ScalarBackend>(0.0,1.0);
  }

  elapsedTime = timer.Elapsed();
  result = sum;

  return elapsedTime;
}

// Vector

double VectorMRG32k3aExp(int nsample, double& result)
{
  // Vector MRG32k3a
  using Double_v = typename VectorBackend::Double_v;
  int vsize = VectorSize<Double_v>();

  vecRng::cxx::MRG32k3a<VectorBackend> rng;
  rng.Initialize();

  static Timer<nanoseconds> timer;
  double elapsedTime = 0.;

  Double_v sum = 0.;

  timer.Start();

  for (int i = 0; i < nsample/vsize ; ++i) {
    sum += rng.Exp<VectorBackend>(1.0);
  }

  elapsedTime = timer.Elapsed();
  for (int i = 0; i < vsize ; ++i) result += sum[i];

  return elapsedTime;
}

double VectorMRG32k3aNormal(int nsample, double& result)
{
  // Vector MRG32k3a
  using Double_v = typename VectorBackend::Double_v;
  int vsize = VectorSize<Double_v>();

  vecRng::cxx::MRG32k3a<VectorBackend> rng;
  rng.Initialize();

  static Timer<nanoseconds> timer;
  double elapsedTime = 0.;

  Double_v sum = 0.;

  timer.Start();

  for (int i = 0; i < nsample/vsize ; ++i) {
    sum += rng.Gauss<VectorBackend>(0.0,1.0);
  }

  elapsedTime = timer.Elapsed();
  for (int i = 0; i < vsize ; ++i) result += sum[i];

  return elapsedTime;
}

} // end namespace vecrng
