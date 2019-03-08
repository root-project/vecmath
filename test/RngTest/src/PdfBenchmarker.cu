#include "PdfBenchmarker.h"
#include "PdfBenchmarker_gpu.h"

#include <cuda.h>
#include <curand_kernel.h>

#include "RngTest.h"

namespace vecrng {

void PdfBenchmarker::RunCuda() {
  int nDevice;
  bool cudaEnabled = false;

  cudaGetDeviceCount(&nDevice);
  if (nDevice > 0) {
    cudaDeviceReset();
    cudaEnabled = true;
  } else {
    printf("Waning: No Cuda Capable Device ...\n");
  }

  // cuda event timing
  cudaEvent_t start;
  cudaEvent_t stop;

  cudaEventCreate(&start);
  cudaEventCreate(&stop);

  // set the default number of threads and thread blocks - should be setable
  // theNThreads should be a power of 2 (due to reduction operations on GPU)
  int theNBlocks = 64;
  int theNThreads = 256;

  // 1. MRG32k3a:

  vecRng::MRG32k3a<vecRng::ScalarBackend> *mrg32k2a =
      new vecRng::MRG32k3a<vecRng::ScalarBackend>();
  vecRng::MRG32k3a<vecRng::ScalarBackend>::State_t *statesMRG32k3a_d = 0;
  cudaMalloc((void **)&statesMRG32k3a_d,
             theNBlocks * theNThreads *
                 sizeof(vecRng::MRG32k3a<vecRng::ScalarBackend>::State_t));
  mrg32k2a->Initialize(statesMRG32k3a_d, theNBlocks * theNThreads);

  // 4 curandStateMRG32k3a
  curandStateMRG32k3a *devStatesMRG32k3a = 0;
  cudaMalloc(&devStatesMRG32k3a,
             theNBlocks * theNThreads * sizeof(curandStateMRG32k3a));
  curand_setup_gpu_gauss(devStatesMRG32k3a, time(NULL), theNBlocks,
                         theNThreads);

  // return values for varification
  double *result_h;
  double *result_c;
  double *result_d;

  result_h = (double *)calloc(theNBlocks, sizeof(double));
  result_c = (double *)calloc(theNBlocks, sizeof(double));
  cudaMalloc((void **)&result_d, theNBlocks * sizeof(double));

  float meanEventTime[kNumberPdf + 1];
  float sigmaEventTime[kNumberPdf + 1];

  double rngEvent[kNumberPdf + 1];

  float *trialEventTime = new float[fRepetition];
  ;

  for (int k = 0; k < kNumberPdf + 1; ++k) {

    meanEventTime[k] = 0.0;
    sigmaEventTime[k] = 0.;
    rngEvent[k] = 0.0;

    float elapsedTotalTime = 0.;

    for (unsigned r = 0; r < fRepetition; ++r) {

      trialEventTime[r] = 0.0;
      cudaMemset(result_d, 0, theNBlocks * theNThreads * sizeof(double));

      if (cudaEnabled) {
        cudaEventRecord(start, 0);

        // call CUDA kernel
        if (k == 0) {
          CudaMRG32k3aGauss(statesMRG32k3a_d, result_d, fNSample, theNBlocks,
                            theNThreads);
        }
        if (k == 1) {
          CurandMRG32k3aGauss(devStatesMRG32k3a, result_d, fNSample, theNBlocks,
                              theNThreads);
        }
        cudaEventRecord(stop, 0);
        cudaEventSynchronize(stop);
        cudaEventElapsedTime(&trialEventTime[r], start, stop);

        // copy the result for varification
        cudaMemcpy(result_h, result_d, theNBlocks * sizeof(double),
                   cudaMemcpyDeviceToHost);

        for (int i = 0; i < theNBlocks; ++i) {
          rngEvent[k] += result_h[i];
        }
        elapsedTotalTime += trialEventTime[r]; // ms
      }
    }

    meanEventTime[k] = elapsedTotalTime / fRepetition;
    float variance = 0;
    for (unsigned r = 0; r < fRepetition; ++r) {
      float delta = (trialEventTime[r] - meanEventTime[k]);
      variance += delta * delta;
    }
    sigmaEventTime[k] = sqrt(variance / fRepetition);
  }
  delete trialEventTime;

  for (int k = 0; k < kNumberPdf + 1; ++k) {
    if (k < kNumberPdf) {
      printf(" %s  CudaBackend   Time = %6.4f +- %6.4f  msec Sum = %g\n",
             PdfName[k], meanEventTime[k], sigmaEventTime[k], rngEvent[k]);
    }
    if (k == kNumberPdf) {
      printf(" %s Nvidia   Time = %6.4f +- %6.4f msec Sum = %g\n",
             "CurandMRG32k3aNormal", meanEventTime[k], sigmaEventTime[k],
             rngEvent[k]);
    }
  }

  // clean up: destory cuda event and free memory on device and host
  cudaEventDestroy(start);
  cudaEventDestroy(stop);

  cudaFree(statesMRG32k3a_d);
  cudaFree(devStatesMRG32k3a);

  cudaFree(result_d);
  free(result_h);
  free(result_c);

  delete mrg32k2a;
}

} // namespace vecrng
