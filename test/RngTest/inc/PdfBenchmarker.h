#ifndef PdfBenchmarker_H
#define PdfBenchmarker_H 1

#include "RngTest.h"

namespace vecrng {

enum PdfIndex { kNullRng = -1, kExp, kGauss, kNumberPdf };
static const char *PdfName[kNumberPdf] = {"MRG32k3aExp", "MRG32k3aNormal"};

class PdfBenchmarker {

public:
  PdfBenchmarker();
  ~PdfBenchmarker();

  int RunBenchmark();

  void SetNSample(const int nsample) { fNSample = nsample; }
  void SetRepetition(const unsigned repetition) { fRepetition = repetition; }

private:
  int RunBenchmarkPdf();

  void RunScalar();
  void RunVector();

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
