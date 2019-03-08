#ifndef RngBenchmarker_H
#define RngBenchmarker_H 1

#include "RngTest.h"

namespace vecrng {

enum RngIndex { kNullRng = -1, kMRG32k3a, kThreefry, kPhilox, kNumberRng };
static const char *RngName[kNumberRng] = {"MRG32k3a", "Threefry", "Philox  "};

#ifdef RNGTEST_MKL
enum VslIndex {
  kNullVsl = -1,
  kMRG32K3A,
  kMCG59,
  kMT17739,
  kSFMT17739,
  kGFSR250,
  kSOBOL32,
  kNumberVsl
};
static const char *VslName[kNumberVsl] = {"MRG32K3A", "MCG59   ", "MT19937 ",
                                          "SMT19937", "GFSR250 ", "SOBOL32 "};
#endif

class RngBenchmarker {

public:
  RngBenchmarker();
  ~RngBenchmarker();

  int RunBenchmark();

  void SetNSample(const int nsample) { fNSample = nsample; }
  void SetRepetition(const unsigned repetition) { fRepetition = repetition; }

private:
  int RunBenchmarkRng();

  void RunTest();
  void RunScalar();
  void RunVector();
  void RunNState();

#ifdef RNGTEST_MKL
  void RunMKLVSL();
#endif

#ifdef RNGTEST_CUDA
  void RunCuda();
#endif

private:
  int fNSample;
  unsigned fRepetition;
  int fVerbosity;
};

} // end namespace vecrng

#endif
