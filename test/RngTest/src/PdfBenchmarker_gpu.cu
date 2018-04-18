#include "VecCoreLib/Rng/MRG32k3a.h"

namespace vecrng {
inline namespace cuda {

const int THREADS_PER_BLOCK = 256;

//do reduction on GPU: THREADS_PER_BLOCK should be 2^N
#define __reduction_on_gpu(Input, Index, Result)     \
{                                                    \
  __syncthreads();                                   \
  int i =  blockDim.x/2;                             \
  while (i != 0) {                                   \
    if( Index < i) Input[Index] += Input[Index + i]; \
    __syncthreads();                                 \
    i/= 2;                                           \
  }                                                  \
  if(sid ==0) Result = Input[0];                     \
} 

__global__
void KernelMRG32k3aGauss(vecRng::MRG32k3a<vecRng::ScalarBackend>::State_t* devStates, double *result, int nsample) 
{
  unsigned int tid = threadIdx.x + blockIdx.x * blockDim.x;
  unsigned int sid = threadIdx.x;

  vecRng::MRG32k3a<vecRng::ScalarBackend> rng(0);

  __shared__ double sum[THREADS_PER_BLOCK];
  double tmp = 0;

  while (tid < nsample) {
    tmp += rng.Gauss<vecRng::ScalarBackend>(&devStates[sid],0.0,1.0);
    tid += blockDim.x * gridDim.x;
  }
  sum[sid] = tmp;

  //do reduction on GPU
  __reduction_on_gpu(sum, sid, result[blockIdx.x]); 
}

//-----------------------------------------------------------------------------
//  Curand MRG32k3a
//-----------------------------------------------------------------------------
__global__
void KernelCurandMRG32k3aGauss(curandStateMRG32k3a *devStates, double *result, int nsample)
{
  unsigned int tid = threadIdx.x + blockIdx.x * blockDim.x;
  unsigned int sid = threadIdx.x;

  curandStateMRG32k3a localState = devStates[sid];

  __shared__ double sum[THREADS_PER_BLOCK];
  double tmp = 0;

  while (tid < nsample) {
    tmp += curand_normal_double(&localState);
    tid += blockDim.x * gridDim.x;
  }
  devStates[sid] = localState;
  sum[sid] = tmp;

  //do reduction on GPU
  __reduction_on_gpu(sum, sid, result[blockIdx.x]); 
}

__global__
void curand_setup_kernel_gauss(curandStateMRG32k3a *devStates, unsigned long seed) 
{
  unsigned int tid = threadIdx.x + blockIdx.x * blockDim.x;
  curand_init(seed, tid, 0, &devStates[tid]);
}


} // end namespace cuda

// Cuda wrapper

  void CudaMRG32k3aGauss(vecRng::MRG32k3a<vecRng::ScalarBackend>::State_t *devStates,
		  double *result,
 		  int nsample,
                  int blocksPerGrid, 
		  int threadsPerBlock) 
{
  KernelMRG32k3aGauss<<<blocksPerGrid, threadsPerBlock>>>(devStates,result,nsample);
}

//-----------------------------------------------------------------------------
//  cuda wrapper for Curand MRG32k3a Kernel
//-----------------------------------------------------------------------------
void CurandMRG32k3aGauss(curandStateMRG32k3a *devStates, 
		         double *result,
 		         int nsample,
                         int blocksPerGrid, 
		         int threadsPerBlock) 
{
  int kstatus = 0;

  KernelCurandMRG32k3aGauss<<<blocksPerGrid, threadsPerBlock>>>(devStates,result,nsample);

  kstatus = cudaThreadSynchronize();
  if(kstatus) std::cout << "MRG32k3a_gpu status = " << kstatus << "\n";
}

void curand_setup_gpu_gauss(curandStateMRG32k3a *devStates, unsigned long seed,  int NBLOCKS, int NTHREADS) {

  int kstatus = 0;

  int threadsPerBlock = NTHREADS;
  int blocksPerGrid = NBLOCKS;

  curand_setup_kernel_gauss<<< blocksPerGrid, threadsPerBlock >>> (devStates,seed);

  kstatus = cudaThreadSynchronize();
  if (kstatus) std::cout << "MRG32k3a: cuda_setup_kernel status = " << kstatus << "\n";
}

} // end namespace vecrng
