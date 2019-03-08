/**
 * RngExample1
 *
 * An example to show how to exchange random states between scalar and vector
 * backend using MRG32k3a
 *
 */

#include <iostream>

#include "RngTest.h"
#include "VecMath/Rng/MRG.h"
#include "VecMath/Rng/MRG32k3a.h"

using namespace vecrng;

void CheckSum(double theory, double sigma, double sample) {
  bool pass = fabs(theory - sample) < sigma;
  if (pass)
    std::cout << " CheckSum PASSED: sampleSum = " << sample << std::endl;
  else
    std::cout << " CheckSum FAILED: sampleSum = " << sample << std::endl;
}

int main(int argc, char *argv[]) {
  // default run with nsample
  int nsample = 100000;
  if (argc >= 2)
    nsample = atoi(argv[1]);
  std::cout << " Number of samples = " << nsample << std::endl;

  // vector size for double pricision
  using Double_v = typename VectorBackend::Double_v;
  int vsize = VectorSize<Double_v>();
  std::cout << " VectorSize<Double_v> = " << vsize << std::endl;

  // using MRG32k3a generator
  typedef vecRng::cxx::MRG32k3a<ScalarBackend> scalarRNG;
  typedef vecRng::cxx::MRG32k3a<VectorBackend> vectorRNG;

  // 1. create a set of scalarRNGs: the size of the given SIMD length
  scalarRNG *srng = new scalarRNG[vsize];
  for (int i = 0; i < vsize; ++i)
    srng[i].Initialize();

  // generate scalar random numbers: total of nsample x vsize
  double ssum = 0;
  for (int i = 0; i < nsample; ++i) {
    for (int i = 0; i < vsize; ++i)
      ssum += srng[i].Uniform<ScalarBackend>();
  }

  // expected mean and sigma of the i.i.d distribution
  double mean = nsample * vsize * 0.5;
  double sigma = sqrt(nsample * vsize / 12.);
  std::cout << " Sample i.i.d [mean,sigma] = [" << mean << "," << sigma << "]"
            << std::endl;

  CheckSum(mean, sigma, ssum);

  // 2. create a vector generator: MRG32k3a<VectorBackend>
  vectorRNG *vrng = new vectorRNG;
  vrng->Initialize();

  // generate vector random numbers
  Double_v vsum(0.);
  for (int i = 0; i < nsample; ++i) {
    vsum += vrng->Uniform<VectorBackend>();
  }

  double avsum = 0.;
  for (int i = 0; i < vsize; ++i) {
    avsum += vsum[i];
  }
  CheckSum(mean, sigma, avsum);

  // 3. gather scalar states to a vector state
  for (int j = 0; j < vsize; ++j) {
    vrng->SetStateAt(j, srng[j].GetState());
  }

  // generate vector random numbers
  Double_v gsum = 0;
  for (int i = 0; i < nsample; ++i) {
    gsum += vrng->Uniform<VectorBackend>();
  }

  double bvsum = 0.;
  for (int i = 0; i < vsize; ++i) {
    bvsum += gsum[i];
  }
  CheckSum(mean, sigma, bvsum);

  // 4. scatter the new vector state back to scalar rngs
  for (int j = 0; j < vsize; ++j) {
    scalarRNG::State_t sstate;
    sstate = vrng->GetStateAt(j);
    srng[j].SetState(&sstate);
  }

  // generate scalar random numbers
  ssum = 0;
  for (int i = 0; i < nsample; ++i) {
    for (int j = 0; j < vsize; ++j)
      ssum += srng[j].Uniform<ScalarBackend>();
  }
  CheckSum(mean, sigma, ssum);

  return 0;
}
