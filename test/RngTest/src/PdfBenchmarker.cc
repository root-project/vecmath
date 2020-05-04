#include "PdfBenchmarker.h"
#include "PdfBenchmarker_cpu.h"

#include "Timer.h"

namespace vecrng {

PdfBenchmarker::PdfBenchmarker() : fNSample(100), fRepetition(1) {}

PdfBenchmarker::~PdfBenchmarker() {}

int PdfBenchmarker::RunBenchmark() {
  printf(" Run PdfBenchmarker with arguments: %d %d\n", fNSample, fRepetition);
  printf(" NSample              = %d\n", fNSample);
  printf(" NRepetitions         = %d\n", fRepetition);
  printf(" VectorSize<Double_v> = %lu\n", VectorSize<Double_v>());

  int errorcode = 0;
  errorcode += RunBenchmarkPdf();
  return (errorcode) ? 1 : 0;
}

int PdfBenchmarker::RunBenchmarkPdf() {
  if (fVerbosity > 0) {
  }

  RunScalar();
  RunVector();

#ifdef RNGTEST_CUDA
  RunCuda();
#endif

  return 0;
}

// Scalar
void PdfBenchmarker::RunScalar() {
  double meanTime[kNumberPdf];
  double sigmaTime[kNumberPdf];
  double resultTotal[kNumberPdf];

  double result = 0;
  double *trialTime = new double[fRepetition];
  ;

  for (unsigned int k = 0; k < kNumberPdf; ++k) {

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

  for (int k = 0; k < kNumberPdf; ++k) {
    printf(" %s  ScalarBackend Time = %4.3f +- %4.3f msec Sum = %g\n",
           PdfName[k], meanTime[k] * 1E-6, sigmaTime[k] * 1E-6, resultTotal[k]);
  }

  delete[] trialTime;
}

// Vector
void PdfBenchmarker::RunVector() {
  double meanTime[kNumberPdf];
  double sigmaTime[kNumberPdf];
  double resultTotal[kNumberPdf];

  double result = 0;
  double *trialTime = new double[fRepetition];
  ;

  for (unsigned int k = 0; k < kNumberPdf; ++k) {

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

  for (int k = 0; k < kNumberPdf; ++k) {
    printf(" %s  VectorBackend Time = %4.3f +- %4.3f msec Sum = %g\n",
           PdfName[k], meanTime[k] * 1E-6, sigmaTime[k] * 1E-6, resultTotal[k]);
  }
}

} // end namespace vecrng
