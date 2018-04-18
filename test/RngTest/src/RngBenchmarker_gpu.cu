#include "VecCoreLib/Rng/MRG32k3a.h"
#include "VecCoreLib/Rng/Threefry.h"
#include "VecCoreLib/Rng/Philox.h"

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
void KernelMRG32k3a(vecRng::MRG32k3a<vecRng::ScalarBackend>::State_t* devStates, double *result, int nsample) 
{
  unsigned int tid = threadIdx.x + blockIdx.x * blockDim.x;
  unsigned int sid = threadIdx.x;

  vecRng::MRG32k3a<vecRng::ScalarBackend> rng(0);

  __shared__ double sum[THREADS_PER_BLOCK];
  double tmp = 0;

  while (tid < nsample) {
    tmp += rng.Uniform<vecRng::ScalarBackend>(&devStates[sid]);
    tid += blockDim.x * gridDim.x;
  }
  sum[sid] = tmp;

  //do reduction on GPU
  __reduction_on_gpu(sum, sid, result[blockIdx.x]); 
}

__global__
void KernelThreefry(vecRng::Threefry<vecRng::ScalarBackend>::State_t* devStates, double *result, int nsample) 
{
  unsigned int tid = threadIdx.x + blockIdx.x * blockDim.x;
  unsigned int sid = threadIdx.x;

  vecRng::Threefry<vecRng::ScalarBackend> rng(0);
 
  __shared__ double sum[THREADS_PER_BLOCK];
  double tmp = 0;

  while (tid < nsample) {
    tmp += rng.Uniform<vecRng::ScalarBackend>(&devStates[sid]);
    tid += blockDim.x * gridDim.x;
  }
  sum[sid] = tmp;

  //do reduction on GPU
  __reduction_on_gpu(sum, sid, result[blockIdx.x]); 
}

__global__
void KernelPhilox(vecRng::Philox<vecRng::ScalarBackend>::State_t* devStates, double *result, int nsample) 
{
  unsigned int tid = threadIdx.x + blockIdx.x * blockDim.x;
  unsigned int sid = threadIdx.x;

  vecRng::Philox<vecRng::ScalarBackend> rng(0);

  __shared__ double sum[THREADS_PER_BLOCK];
  double tmp = 0;

  while (tid < nsample) {
    tmp += rng.Uniform<vecRng::ScalarBackend>(&devStates[sid]);
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
void KernelCurandMRG32k3a(curandStateMRG32k3a *devStates, double *result, int nsample)
{
  unsigned int tid = threadIdx.x + blockIdx.x * blockDim.x;
  unsigned int sid = threadIdx.x;

  curandStateMRG32k3a localState = devStates[sid];

  __shared__ double sum[THREADS_PER_BLOCK];
  double tmp = 0;

  while (tid < nsample) {
    tmp += curand_uniform_double(&localState);
    tid += blockDim.x * gridDim.x;
  }
  devStates[sid] = localState;
  sum[sid] = tmp;

  //do reduction on GPU
  __reduction_on_gpu(sum, sid, result[blockIdx.x]); 
}

__global__
void curand_setup_kernel(curandStateMRG32k3a *devStates, unsigned long seed) 
{
  unsigned int tid = threadIdx.x + blockIdx.x * blockDim.x;
  curand_init(seed, tid, 0, &devStates[tid]);
}

//-----------------------------------------------------------------------------
//  Curand Philox
//-----------------------------------------------------------------------------
__global__
void KernelCurandPhilox(curandStatePhilox4_32_10_t *devStates, double *result, int nsample)
{
  unsigned int tid = threadIdx.x + blockIdx.x * blockDim.x;
  unsigned int sid = threadIdx.x;

  curandStatePhilox4_32_10_t localState = devStates[tid];

  __shared__ double sum[THREADS_PER_BLOCK];
  double tmp = 0;

  while (tid < nsample) {
    tmp += curand_uniform_double(&localState);
    tid += blockDim.x * gridDim.x;
  }
  devStates[sid] = localState;
  sum[sid] = tmp;

  //do reduction on GPU
  __reduction_on_gpu(sum, sid, result[blockIdx.x]); 
}

__global__
void curand_setup_kernel(curandStatePhilox4_32_10_t *devStates, unsigned long seed) 
{
  unsigned int tid = threadIdx.x + blockIdx.x * blockDim.x;
  curand_init(seed, tid, 0, &devStates[tid]);
}


} // end namespace cuda

// Cuda wrapper

void CudaMRG32k3a(vecRng::MRG32k3a<vecRng::ScalarBackend>::State_t *devStates,
		  double *result,
 		  int nsample,
                  int blocksPerGrid, 
		  int threadsPerBlock) 
{
  KernelMRG32k3a<<<blocksPerGrid, threadsPerBlock>>>(devStates,result,nsample);
}

void CudaThreefry(vecRng::Threefry<vecRng::ScalarBackend>::State_t *devStates,
		  double *result,
 		  int nsample,
                  int blocksPerGrid, 
		  int threadsPerBlock) 
{
   KernelThreefry<<<blocksPerGrid, threadsPerBlock>>>(devStates,result,nsample);
}

void CudaPhilox(vecRng::Philox<vecRng::ScalarBackend>::State_t *devStates,
		double *result,
 		int nsample,
                int blocksPerGrid, 
		int threadsPerBlock) 
{
   KernelPhilox<<<blocksPerGrid, threadsPerBlock>>>(devStates,result,nsample);
}

//-----------------------------------------------------------------------------
//  cuda wrapper for Curand MRG32k3a Kernel
//-----------------------------------------------------------------------------
void CurandMRG32k3a(curandStateMRG32k3a *devStates, 
		  double *result,
 		  int nsample,
                  int blocksPerGrid, 
		  int threadsPerBlock) 
{
  int kstatus = 0;

  KernelCurandMRG32k3a<<<blocksPerGrid, threadsPerBlock>>>(devStates,result,nsample);

  kstatus = cudaThreadSynchronize();
  if(kstatus) std::cout << "MRG32k3a_gpu status = " << kstatus << "\n";
}

void curand_setup_gpu(curandStateMRG32k3a *devStates, unsigned long seed,  int NBLOCKS, int NTHREADS) {

  int kstatus = 0;

  int threadsPerBlock = NTHREADS;
  int blocksPerGrid = NBLOCKS;

  curand_setup_kernel<<< blocksPerGrid, threadsPerBlock >>> (devStates,seed);

  kstatus = cudaThreadSynchronize();
  if (kstatus) std::cout << "MRG32k3a: cuda_setup_kernel status = " << kstatus << "\n";
}

//-----------------------------------------------------------------------------
//  cuda wrapper for Curand Philox4_32_10 Kernel
//-----------------------------------------------------------------------------
void CurandPhilox(curandStatePhilox4_32_10_t *devStates, 
		  double *result,
 		  int nsample,
                  int blocksPerGrid, 
		  int threadsPerBlock) 
{
  int kstatus = 0;

  KernelCurandPhilox<<<blocksPerGrid, threadsPerBlock>>>(devStates,result,nsample);

  kstatus = cudaThreadSynchronize();
  if(kstatus) std::cout << "CurandPhilox status = " << kstatus << "\n";
}

void curand_setup_gpu(curandStatePhilox4_32_10_t *devStates, unsigned long seed,  int NBLOCKS, int NTHREADS) {

  int kstatus = 0;

  int threadsPerBlock = NTHREADS;
  int blocksPerGrid = NBLOCKS;

  curand_setup_kernel<<< blocksPerGrid, threadsPerBlock >>> (devStates,seed);

  kstatus = cudaThreadSynchronize();
  if (kstatus) std::cout << "Philox: cuda_setup_kernel status = " << kstatus << "\n";
}

} // end namespace vecrng
