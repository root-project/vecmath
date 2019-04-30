#include "RngBenchmarker.h"
#include "RngBenchmarker_cpu.h"

#include "Timer.h"

namespace vecrng {

RngBenchmarker::RngBenchmarker() : fNSample(100), fRepetition(1) {}

RngBenchmarker::~RngBenchmarker() {}

int RngBenchmarker::RunBenchmark() {
  printf(" Run RngBenchmarker with arguments: %d %d\n", fNSample, fRepetition);
  printf(" NSample              = %d\n", fNSample);
  printf(" NRepetitions         = %d\n", fRepetition);
  printf(" VectorSize<UInt32_v> = %lu\n", VectorSize<UInt32_v>());
  printf(" VectorSize<UInt64_v> = %lu\n", VectorSize<UInt64_v>());
  printf(" VectorSize<Double_v> = %lu\n", VectorSize<Double_v>());

  int errorcode = 0;
  errorcode += RunBenchmarkRng();
  return (errorcode) ? 1 : 0;
}

int RngBenchmarker::RunBenchmarkRng() {
  if (fVerbosity > 0) {
  }

  RunTest();
  RunScalar();
  RunVector();
  RunNState();

#ifdef RNGTEST_MKL
  RunMKLVSL();
#endif

#ifdef RNGTEST_CUDA
  RunCuda();
#endif

  return 0;
}

// Reference
void RngBenchmarker::RunTest() {
  // test GNU rand() as a reference

  double elapsedTotal = 0.;
  double resultTotal = 0.;

  double *trialTime = new double[fRepetition];
  ;

  static Timer<nanoseconds> timer;

  for (unsigned r = 0; r < fRepetition; ++r) {

    double result = 0.;

    timer.Start();

    for (int i = 0; i < fNSample; ++i) {
      result += (double)rand() / RAND_MAX;
    }

    trialTime[r] = timer.Elapsed();

    elapsedTotal += trialTime[r];
    resultTotal += result;
  }

  double meanTime = elapsedTotal / fRepetition;
  double variance = 0;

  for (unsigned r = 0; r < fRepetition; ++r) {
    double delta = (trialTime[r] - meanTime);
    variance += delta * delta;
  }
  double sigmaTime = sqrt(variance / fRepetition);

  delete[] trialTime;

  printf(" %s std::rand()   Time = %4.3f +- %4.3f msec Sum = %g\n", "TestRand ",
         meanTime * 1E-6, sigmaTime * 1E-6, resultTotal);
}

// Scalar
void RngBenchmarker::RunScalar() {
  double meanTime[kNumberRng];
  double sigmaTime[kNumberRng];
  double resultTotal[kNumberRng];

  double result = 0;
  double *trialTime = new double[fRepetition];
  ;

  for (unsigned int k = 0; k < kNumberRng; ++k) {

    meanTime[k] = 0.;
    sigmaTime[k] = 0.;
    resultTotal[k] = 0.;

    double elapsedTotal = 0.0;

    for (unsigned r = 0; r < fRepetition; ++r) {
      trialTime[r] = 0.0;
      result = 0.0;
      trialTime[r] = ScalarKernelFunc[k](fNSample, result);
      elapsedTotal += trialTime[r];
      resultTotal[k] += result;
    }

    meanTime[k] = elapsedTotal / fRepetition;
    double variance = 0;

    for (unsigned r = 0; r < fRepetition; ++r) {
      double delta = (trialTime[r] - meanTime[k]);
      variance += delta * delta;
    }
    sigmaTime[k] = sqrt(variance / fRepetition);
  }

  for (int k = 0; k < kNumberRng; ++k) {
    printf(" %s  ScalarBackend Time = %4.3f +- %4.3f msec Sum = %g\n",
           RngName[k], meanTime[k] * 1E-6, sigmaTime[k] * 1E-6, resultTotal[k]);
  }

  delete[] trialTime;
}

// Vector
void RngBenchmarker::RunVector() {
  double meanTime[kNumberRng];
  double sigmaTime[kNumberRng];
  double resultTotal[kNumberRng];

  double result = 0;
  double *trialTime = new double[fRepetition];
  ;

  for (unsigned int k = 0; k < kNumberRng; ++k) {

    meanTime[k] = 0.;
    sigmaTime[k] = 0.;
    resultTotal[k] = 0.;

    double elapsedTotal = 0.0;

    for (unsigned r = 0; r < fRepetition; ++r) {
      trialTime[r] = 0.0;
      result = 0.0;
      trialTime[r] = VectorKernelFunc[k](fNSample, result);
      elapsedTotal += trialTime[r];
      resultTotal[k] += result;
    }

    meanTime[k] = elapsedTotal / fRepetition;
    double variance = 0;

    for (unsigned r = 0; r < fRepetition; ++r) {
      double delta = (trialTime[r] - meanTime[k]);
      variance += delta * delta;
    }
    sigmaTime[k] = sqrt(variance / fRepetition);
  }
  delete[] trialTime;

  for (int k = 0; k < kNumberRng; ++k) {
    printf(" %s  VectorBackend Time = %4.3f +- %4.3f msec Sum = %g\n",
           RngName[k], meanTime[k] * 1E-6, sigmaTime[k] * 1E-6, resultTotal[k]);
  }
}

void RngBenchmarker::RunNState() {
  double meanTime[kNumberRng];
  double sigmaTime[kNumberRng];
  double resultTotal[kNumberRng];

  double result = 0;
  double *trialTime = new double[fRepetition];
  ;

  for (unsigned int k = 0; k < kNumberRng; ++k) {

    meanTime[k] = 0.;
    sigmaTime[k] = 0.;
    resultTotal[k] = 0.;

    double elapsedTotal = 0.0;

    for (unsigned r = 0; r < fRepetition; ++r) {
      trialTime[r] = 0.0;
      result = 0.0;
      trialTime[r] = StateKernelFunc[k](fNSample, result);
      elapsedTotal += trialTime[r];
      resultTotal[k] += result;
    }

    meanTime[k] = elapsedTotal / fRepetition;
    double variance = 0;

    for (unsigned r = 0; r < fRepetition; ++r) {
      double delta = (trialTime[r] - meanTime[k]);
      variance += delta * delta;
    }
    sigmaTime[k] = sqrt(variance / fRepetition);
  }
  delete[] trialTime;

  for (int k = 0; k < kNumberRng; ++k) {
    printf(" %s  ScalarNstates Time = %4.3f +- %4.3f msec Sum = %g\n",
           RngName[k], meanTime[k] * 1E-6, sigmaTime[k] * 1E-6, resultTotal[k]);
  }
}

#ifdef RNGTEST_MKL
void RngBenchmarker::RunMKLVSL() {
  double meanTime[kNumberVsl];
  double sigmaTime[kNumberVsl];
  double resultTotal[kNumberVsl];

  double result = 0;
  double *trialTime = new double[fRepetition];
  ;

  for (unsigned int k = 0; k < kNumberVsl; ++k) {

    meanTime[k] = 0.;
    sigmaTime[k] = 0.;
    resultTotal[k] = 0.;

    double elapsedTotal = 0.0;

    for (unsigned r = 0; r < fRepetition; ++r) {
      trialTime[r] = 0.0;
      result = 0.0;
      trialTime[r] = VSLKernelFunc[k](fNSample, result);
      elapsedTotal += trialTime[r];
      resultTotal[k] += result;
    }

    meanTime[k] = elapsedTotal / fRepetition;
    double variance = 0;

    for (unsigned r = 0; r < fRepetition; ++r) {
      double delta = (trialTime[r] - meanTime[k]);
      variance += delta * delta;
    }
    sigmaTime[k] = sqrt(variance / fRepetition);
  }
  delete[] trialTime;

  for (int k = 0; k < kNumberVsl; ++k) {
    printf(" %s  Intel-MKL/VSL Time = %3.3f +- %3.3f msec Sum = %g\n",
           VslName[k], meanTime[k] * 1E-6, sigmaTime[k] * 1E-6, resultTotal[k]);
  }
}
#endif

} // end namespace vecrng
