#include "RngBenchmarker.h"
#include "RngBenchmarker_gpu.h"

#include <cuda.h>
#include <curand_kernel.h>

#include "RngTest.h"

namespace vecrng {

void RngBenchmarker::RunCuda()
{
  int nDevice;
  bool cudaEnabled = false;

  cudaGetDeviceCount(&nDevice);
  if(nDevice > 0) {
    cudaDeviceReset();
    cudaEnabled = true;
  }
  else {
    printf("Waning: No Cuda Capable Device ...\n");
  }

  //cuda event timing
  cudaEvent_t start;
  cudaEvent_t stop;

  cudaEventCreate (&start);
  cudaEventCreate (&stop);

  //set the default number of threads and thread blocks - should be setable
  //theNThreads should be a power of 2 (due to reduction operations on GPU)
  int theNBlocks  =   64;
  int theNThreads =  256;

  //1. MRG32k3a:

  vecRng::MRG32k3a<vecRng::ScalarBackend> *mrg32k2a = new vecRng::MRG32k3a<vecRng::ScalarBackend>();
  vecRng::MRG32k3a<vecRng::ScalarBackend>::State_t* statesMRG32k3a_d = 0; 
  cudaMalloc((void**)&statesMRG32k3a_d, theNBlocks*theNThreads*sizeof(vecRng::MRG32k3a<vecRng::ScalarBackend>::State_t));
  mrg32k2a->Initialize(statesMRG32k3a_d, theNBlocks*theNThreads);

  //2. Threefry:
  vecRng::Threefry<vecRng::ScalarBackend> *threefry = new vecRng::Threefry<vecRng::ScalarBackend>();
  vecRng::Threefry<vecRng::ScalarBackend>::State_t* statesThreefry_d = 0; 
  cudaMalloc((void**)&statesThreefry_d, theNBlocks*theNThreads*sizeof(vecRng::Threefry<vecRng::ScalarBackend>::State_t));
  threefry->Initialize(statesThreefry_d, theNBlocks*theNThreads);

  //Philox:
  vecRng::Philox<vecRng::ScalarBackend> *philox = new vecRng::Philox<vecRng::ScalarBackend>();
  vecRng::Philox<vecRng::ScalarBackend>::State_t* statesPhilox_d = 0; 
  cudaMalloc((void**)&statesPhilox_d, theNBlocks*theNThreads*sizeof(vecRng::Philox<vecRng::ScalarBackend>::State_t));
  philox->Initialize(statesPhilox_d, theNBlocks*theNThreads);

  //4 curandStateMRG32k3a
  curandStateMRG32k3a* devStatesMRG32k3a = 0;
  cudaMalloc(&devStatesMRG32k3a,theNBlocks*theNThreads*sizeof(curandStateMRG32k3a));
  curand_setup_gpu(devStatesMRG32k3a, time(NULL), theNBlocks, theNThreads);

  //4 curandStatePhilox
  curandStatePhilox4_32_10_t* devStatesPhilox = 0;
  cudaMalloc(&devStatesPhilox,theNBlocks*theNThreads*sizeof(curandStateMRG32k3a));
  curand_setup_gpu(devStatesPhilox, time(NULL), theNBlocks, theNThreads);

  //return values for varification
  double *result_h;
  double *result_c;
  double *result_d;

  result_h = (double*) calloc(theNBlocks, sizeof(double));
  result_c = (double*) calloc(theNBlocks, sizeof(double));
  cudaMalloc((void**)&result_d,theNBlocks*sizeof(double));

  float meanEventTime[kNumberRng +2];
  float sigmaEventTime[kNumberRng +2];

  double rngEvent[kNumberRng +2];

  float *trialEventTime = new float [fRepetition];;

  for (int k = 0; k < kNumberRng + 2; ++k) {

    meanEventTime[k] = 0.0;
    sigmaEventTime[k] = 0.;
    rngEvent[k] = 0.0;

    float elapsedTotalTime = 0.;

    for (unsigned r = 0; r < fRepetition; ++r) {

      trialEventTime[r] = 0.0;
      cudaMemset(result_d,0,theNBlocks*theNThreads*sizeof(double));
  
      if(cudaEnabled) {
        cudaEventRecord (start,0);

        //call CUDA kernel
	if(k == 0) {
	  CudaMRG32k3a(statesMRG32k3a_d, result_d, fNSample, theNBlocks, theNThreads);
        }
	if(k == 1) {
	  CudaThreefry(statesThreefry_d, result_d, fNSample, theNBlocks, theNThreads);
	}
	if(k == 2) {
	  CudaPhilox(statesPhilox_d, result_d, fNSample, theNBlocks, theNThreads);
	}

	if(k == 3) {
          CurandMRG32k3a(devStatesMRG32k3a,result_d,fNSample,theNBlocks,theNThreads);
	}

	if(k == 4) {
          CurandPhilox(devStatesPhilox,result_d,fNSample,theNBlocks,theNThreads);
	}

        cudaEventRecord (stop,0);
        cudaEventSynchronize (stop);
        cudaEventElapsedTime (&trialEventTime[r],start,stop);

        //copy the result for varification
        cudaMemcpy(result_h,result_d,theNBlocks*sizeof(double),cudaMemcpyDeviceToHost);
         
        for(int i = 0 ; i < theNBlocks ; ++i) {
          rngEvent[k] += result_h[i]; 
        }
        elapsedTotalTime += trialEventTime[r]; //ms
      }
    }

    meanEventTime[k] = elapsedTotalTime/fRepetition;
    float variance = 0;
    for (unsigned r = 0; r < fRepetition; ++r) {
      float delta  = (trialEventTime[r] - meanEventTime[k]);
      variance += delta*delta;
    }
    sigmaEventTime[k] = sqrt(variance/fRepetition);   
  }
  delete trialEventTime;

  for (int k = 0; k < kNumberRng + 2; ++k) {
    if(k < kNumberRng) {
      printf(" %s  CudaBackend   Time = %6.4f +- %6.4f  msec Sum = %g\n", 
   	     RngName[k], meanEventTime[k], sigmaEventTime[k], rngEvent[k]);
    }
    if(k== kNumberRng) {
      printf(" %s Nvidia   Time = %6.4f +- %6.4f msec Sum = %g\n", 
   	     "CurandMRG32k3a", meanEventTime[k], sigmaEventTime[k], rngEvent[k]);
    }
    if(k== kNumberRng+1) {
      printf(" %s Nvidia   Time = %6.4f +- %6.4f msec Sum = %g\n", 
   	     "CurandPhilox  ", meanEventTime[k], sigmaEventTime[k], rngEvent[k]);
    }
  }

  //clean up: destory cuda event and free memory on device and host
  cudaEventDestroy(start);
  cudaEventDestroy(stop);

  cudaFree(statesMRG32k3a_d);
  cudaFree(statesThreefry_d);
  cudaFree(statesPhilox_d);
  cudaFree(devStatesMRG32k3a);
  cudaFree(devStatesPhilox);

  cudaFree(result_d);
  free(result_h);
  free(result_c);

  delete mrg32k2a;
  delete threefry;
  delete philox;
}

} // end of vecrng namespace

