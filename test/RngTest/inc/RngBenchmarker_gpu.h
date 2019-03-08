#ifndef RNG_BENCHMARKER_GPU_H
#define RNG_BENCHMARKER_GPU_H 1

#include "VecMath/Rng/MRG32k3a.h"
#include "VecMath/Rng/Philox.h"
#include "VecMath/Rng/Threefry.h"

#include <curand_kernel.h>

namespace vecrng {

// Cuda Wrapper
void CudaMRG32k3a(vecRng::MRG32k3a<vecRng::ScalarBackend>::State_t *devStates,
                  double *result, int nsample, int blocksPerGrid,
                  int threadsPerBlock);

void CudaThreefry(vecRng::Threefry<vecRng::ScalarBackend>::State_t *devStates,
                  double *result, int nsample, int blocksPerGrid,
                  int threadsPerBlock);

void CudaPhilox(vecRng::Philox<vecRng::ScalarBackend>::State_t *devStates,
                double *result, int nsample, int blocksPerGrid,
                int threadsPerBlock);

// Curand Test
void CurandMRG32k3a(curandStateMRG32k3a *devStates, double *result, int nsample,
                    int blocksPerGrid, int threadsPerBlock);

void curand_setup_gpu(curandStateMRG32k3a *devStates, unsigned long seed,
                      int blocksPerGrid, int threadsPerBlock);

void CurandPhilox(curandStatePhilox4_32_10_t *devStates, double *result,
                  int nsample, int blocksPerGrid, int threadsPerBlock);

void curand_setup_gpu(curandStatePhilox4_32_10_t *devStates, unsigned long seed,
                      int blocksPerGrid, int threadsPerBlock);

} // end namespace vecrng

#endif
