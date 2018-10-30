#ifndef VECRNG_PHILOX_H
#define VECRNG_PHILOX_H 1

/**
  Philox: A SIMD/SIMT implementation of Philox (4x32-10) based on philox.h 
          of Random123 (version 1.09)
  
  Copyright 2010-2011, D. E. Shaw Research.
  All rights reserved.

  Redistribution and use in source and binary forms, with or without
  modification, are permitted provided that the following conditions are
  met:
  
  * Redistributions of source code must retain the above copyright
    notice, this list of conditions, and the following disclaimer.
  
  * Redistributions in binary form must reproduce the above copyright
    notice, this list of conditions, and the following disclaimer in the
    documentation and/or other materials provided with the distribution.
  
  * Neither the name of D. E. Shaw Research nor the names of its
    contributors may be used to endorse or promote products derived from
    this software without specific prior written permission.
  
  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
  "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
  LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
  A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
  OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
  SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
  LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
  DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
  THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
  (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
  OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

#include "R123.h"
#include "VecRNG.h"

#include <limits.h>

namespace vecRng {
inline namespace VECRNG_IMPL_NAMESPACE {

template <typename BackendT> class Philox;

// struct Philox (random state of Philox-4x32-10)

template <typename BackendT>
struct RNG_traits<Philox<BackendT> > {
  struct State {
    R123::array_t<BackendT,4> ctr;
    R123::array_t<BackendT,2> key;
    R123::array_t<BackendT,4> ukey;
    unsigned int index;
  };
  using State_t = State;
}; 

template <typename BackendT>
class Philox : public VecRNG<Philox<BackendT> > {

public:
  using State_t = typename RNG_traits<Philox<BackendT> >::State_t;
  using State_s = typename RNG_traits<Philox<ScalarBackend> >::State_t;
  using State_v = typename RNG_traits<Philox<VectorBackend> >::State_t;

private:
  static long long fSeed;

public:

  VECCORE_ATT_HOST_DEVICE
  Philox() {} 

  VECCORE_ATT_HOST_DEVICE
  Philox(State_t *states) : VecRNG<Philox<BackendT> >(states) {}

  VECCORE_ATT_HOST_DEVICE
  ~Philox() {}

  VECCORE_ATT_HOST_DEVICE
  Philox(const Philox &rng); 

  // Mandatory methods - static inheritance

  VECCORE_ATT_HOST
  VECCORE_FORCE_INLINE
  void Initialize();

  // Initialize with a unique stream number
  VECCORE_ATT_HOST
  VECCORE_FORCE_INLINE
  void Initialize(long streamId);
  
  // Initialize a set of states of which size is nthreads
  VECCORE_ATT_HOST
  VECCORE_FORCE_INLINE
  void Initialize(State_t *states, unsigned int nthreads);
  
  // Returns pRNG<BackendT> between 0 and 1 (excluding the end points).
  template <typename ReturnTypeBackendT>
  VECCORE_ATT_HOST_DEVICE 
  typename ReturnTypeBackendT::Double_v Kernel(State_t& state);

  //Copy a scalar state explicitly to the i-th lane of the vector state
  VECCORE_FORCE_INLINE
  VECCORE_ATT_HOST_DEVICE
  void SetStateAt(unsigned int i, State_s *state); 

  //Return the ith-lane of the vector state to a scalar state
  VECCORE_FORCE_INLINE
  VECCORE_ATT_HOST_DEVICE
  State_s GetStateAt(unsigned int i);

  // Auxiliary methods

  VECCORE_ATT_HOST
  VECCORE_FORCE_INLINE
  void AdvanceState(long long n);
  
  VECCORE_ATT_HOST_DEVICE void SetSeed(long long seed) { fSeed = seed; }

  VECCORE_ATT_HOST_DEVICE long long GetSeed() const { return fSeed; }

  VECCORE_FORCE_INLINE
  VECCORE_ATT_HOST void PrintState() const { return PrintState(*(this->fState)); }

  VECCORE_ATT_HOST void PrintState(State_t const &state) const;

private:
  // the mother is friend of this
  friend class VecRNG<Philox<BackendT> >;

  // Set the stream to the next stream/substream.
  VECCORE_ATT_HOST inline void SetNextStream(State_t *state);

  VECCORE_ATT_HOST inline void SetNextSubstream();

  // Increase counter
  VECCORE_ATT_HOST_DEVICE 
  inline void IncreaseCounter(State_t *state);
  
  // Philox utility methods
  VECCORE_ATT_HOST_DEVICE
  inline  typename BackendT::UInt32_v Mulhilo32(typename BackendT::UInt32_v a, 
                                                typename BackendT::UInt32_v b, 
                                                typename BackendT::UInt32_v *hip);

  VECCORE_ATT_HOST_DEVICE inline void Philox4x32bumpkey(R123::array_t<BackendT,2> key);

  VECCORE_ATT_HOST_DEVICE void Philox4x32round(R123::array_t<BackendT,4> crt, 
  					       R123::array_t<BackendT,2> key);

  VECCORE_ATT_HOST_DEVICE void Gen(R123::array_t<BackendT,4> ctr, R123::array_t<BackendT,2> key,
                                   R123::array_t<BackendT,4> out);

};

// The default seed of Philox
template <class BackendT> long long Philox<BackendT>::fSeed = 12345;

//
// Class Implementation
//  

// Copy constructor
template <typename BackendT>
VECCORE_ATT_HOST_DEVICE
Philox<BackendT>::Philox(const Philox<BackendT> &rng) : VecRNG<Philox<BackendT> >()
{
  this->fState->index = rng.fState->index;
  for(size_t i = 0 ; i < 4 ; ++i) {
    this->fState->ctr[i] = rng.fState->ctr[i];
    if(i < 2) this->fState->key[i] = rng.fState->key[i];
    this->fState->ukey[i]= rng.fState->ukey[i];
  }  
  fSeed = rng.fSeed; 
}

// Set a new set of keys for the vector of next streams using the unique seed
template <typename BackendT>
VECCORE_ATT_HOST inline void Philox<BackendT>::SetNextStream (State_t *state)
{
  for(size_t i = 0 ; i < VectorSize<UInt32_v>() ; ++i) { 
    for(size_t j = 0 ; j < 4 ; ++j) state->ctr[j][i] = 0;
    state->key[0][i] = (unsigned int)(fSeed);
    state->key[1][i] = (unsigned int)(fSeed>>32);
    ++fSeed;
  }
}

// Scalar Specialization of SetNextStream
template <>
VECCORE_ATT_HOST inline void Philox<ScalarBackend>::SetNextStream(State_s *state)
{
  for(size_t i = 0 ; i < 4 ; ++i) state->ctr[i] = 0;
  state->key[0] = (unsigned int)(fSeed);
  state->key[1] = (unsigned int)(fSeed>>32);
  ++fSeed;
}

// Reset the current stream to the next substream - skipahead 
template <typename BackendT>
VECCORE_ATT_HOST inline void Philox<BackendT>::SetNextSubstream ()
{
  for(size_t i = 0 ; i < VectorSize<UInt32_v>() ; ++i) { 
    unsigned int nlo = (unsigned int)(R123::PHILOX_SKIP_AHEAD);
    unsigned int nhi = (unsigned int)(R123::PHILOX_SKIP_AHEAD>>32);

    this->fState->ctr[0] += nlo;
    if( this->fState->ctr[0] < nlo ) nhi++;
    this->fState->ctr[1] += nhi;
    if(nhi <= this->fState->ctr[1]) continue;
    if(++this->fState->ctr[2]) continue;
    ++this->fState->ctr[3];
  }
}

// Scalar specialization of SetNextSubstream
template <>
VECCORE_ATT_HOST inline void Philox<ScalarBackend>::SetNextSubstream ()
{
  unsigned int nlo = (unsigned int)(R123::PHILOX_SKIP_AHEAD);
  unsigned int nhi = (unsigned int)(R123::PHILOX_SKIP_AHEAD>>32);

  this->fState->ctr[0] += nlo;
  if( this->fState->ctr[0] < nlo ) nhi++;
  this->fState->ctr[1] += nhi;
  if(nhi <= this->fState->ctr[1]) return;
  if(++this->fState->ctr[2]) return;
  ++this->fState->ctr[3];
}

// Set the seed for each stream of SIMD to the starting state of each substream
template <class BackendT>
inline VECCORE_ATT_HOST void Philox<BackendT>::Initialize()
{
  //set initial counter and key
  this->fState->index = 0;
  SetNextStream(this->fState);
}

template <typename BackendT>
VECCORE_ATT_HOST
VECCORE_FORCE_INLINE
void Philox<BackendT>::Initialize(long streamId)
{
  this->fState->index = 0;

  //start from the default state and skip to the streamId stream
  fSeed = 12345 + streamId*VectorSize<UInt32_v>();

  for(size_t i = 0 ; i < VectorSize<UInt32_v>() ; ++i) { 
    for(size_t j = 0 ; j < 4 ; ++j) this->fState->ctr[j][i] = 0;
    this->fState->key[0][i] = (unsigned int)(fSeed);
    this->fState->key[1][i] = (unsigned int)(fSeed>>32);
    ++fSeed;
  }
}

template <>
VECCORE_ATT_HOST
VECCORE_FORCE_INLINE
void Philox<ScalarBackend>::Initialize(long streamId)
{
  this->fState->index = 0;
  
  //start from the default state and skip to the streamId-th stream
  fSeed = 12345 + streamId;

  for(size_t j = 0 ; j < 4 ; ++j) this->fState->ctr[j] = 0;
  this->fState->key[0] = (unsigned int)(fSeed);
  this->fState->key[1] = (unsigned int)(fSeed>>32);
  ++fSeed;
}
 
// Specialization of Initialize for SIMT
template <>
VECCORE_ATT_HOST 
VECCORE_FORCE_INLINE
void Philox<ScalarBackend>::Initialize(State_s *states, unsigned int nthreads)
{
  State_s* hstates 
   = (State_s *) malloc (nthreads*sizeof(State_s));

  for (unsigned int tid = 0 ; tid < nthreads ; ++tid) {
    //initialize initial seed/state by the unique tid number
    hstates[tid].index = 0; 
    fSeed = tid;
    SetNextStream(&hstates[tid]);
  }
#ifdef VECCORE_CUDA
  cudaMemcpy(states, hstates, nthreads*sizeof(State_s), cudaMemcpyHostToDevice);
#else
  memcpy(states, hstates, nthreads*sizeof(State_s));
#endif
  free(hstates);
}

// Increase counter of each element of the counter (ctr) vector
template <typename BackendT>
VECCORE_ATT_HOST_DEVICE 
inline void Philox<BackendT>::IncreaseCounter(State_t *state)
{
  size_t vsize = VectorSize<UInt32_v>();
  for(size_t iv = 0 ; iv < vsize ; ++iv) {
    if( ++state->ctr[0][iv]) continue;
    if( ++state->ctr[1][iv]) continue;
    if( ++state->ctr[2][iv]) continue;
    ++state->ctr[3][iv];
  }
}

// Increase counter of each element of the counter (ctr) vector
template <>
VECCORE_ATT_HOST_DEVICE 
inline void Philox<ScalarBackend>::IncreaseCounter(State_s *state)
{
  if( ++state->ctr[0]) return;
  if( ++state->ctr[1]) return;
  if( ++state->ctr[2]) return;
  ++state->ctr[3];
}

// Print information of the current state
template <typename BackendT>
VECCORE_ATT_HOST void Philox<BackendT>::PrintState(State_t const& state) const
{
  std::cout << "index = " << state.index << std::endl;
  for(size_t i = 0 ; i < 2 ; ++i) {
    std::cout << "key[" << i << "] = " <<  state.key[i] << std::endl;
  }
}

template <class BackendT>
VECCORE_FORCE_INLINE
VECCORE_ATT_HOST void Philox<BackendT>::AdvanceState(long long n)
{
  unsigned int nlo = (unsigned int)(n);
  unsigned int nhi = (unsigned int)(n>>32);

  for(size_t i = 0 ; i < VectorSize<UInt32_v>() ; ++i) { 
    this->fState->ctr[0][i] += nlo;
    if(this->fState->ctr[0][i]) nhi++;
    this->fState->ctr[1][i] += nhi;
    if(nhi <= this->fState->ctr[1][i]) return;
    if(++(this->fState->ctr[2][i])) return;
    ++(this->fState->ctr[3][i]);
  }
}

template <>
VECCORE_FORCE_INLINE
VECCORE_ATT_HOST void Philox<ScalarBackend>::AdvanceState(long long n)
{
  unsigned int nlo = (unsigned int)(n);
  unsigned int nhi = (unsigned int)(n>>32);

  this->fState->ctr[0] += nlo;
  if(this->fState->ctr[0]) nhi++;
  this->fState->ctr[1] += nhi;
  if(nhi <= this->fState->ctr[1]) return;
  if(++(this->fState->ctr[2])) return;
  ++(this->fState->ctr[3]);
}

// Kernel to generate a vector(scalar) of next random number(s)
template <class BackendT>
template <class ReturnTypeBackendT>
VECCORE_ATT_HOST_DEVICE  
typename ReturnTypeBackendT::Double_v Philox<BackendT>::Kernel(State_t& state)
{
  using Double_v = typename ReturnTypeBackendT::Double_v;
  Double_v u(0.0);

  if(state.index == 0 ) {
    //get a new state and generate 128 bits of pseudo randomness 
    Gen(state.ctr,state.key,state.ukey);

    //construct 8xUInt32 (ukey) to 4xUInt64, then convert to Double_v using UINT32_MAX
    u = static_cast<Double_v>( (state.ukey[state.index]) )/UINT32_MAX;
    
    //state index and increase counter
    ++(state.index);
    IncreaseCounter(&state);
  }
  else {  
    //    u = (Double_v)(state.ukey[state.index] >>11 ) * a + b;
    u = static_cast<Double_v>( (state.ukey[state.index]) )/UINT32_MAX; 
    ++state.index;
    if(state.index == 4) state.index = 0;
  }
  
  return u;
}

// Sepecialization for the scalar method of VectorBackend pRNG class
#ifndef VECCORE_CUDA
template <>
template <>
VECCORE_FORCE_INLINE
VECCORE_ATT_HOST_DEVICE ScalarBackend::Double_v 
Philox<VectorBackend>::Kernel<ScalarBackend>(State_v& state)
{
  return Kernel<VectorBackend>(state)[0]; 
}
#endif

// Philox utility methods
template <class BackendT>
inline
VECCORE_ATT_HOST_DEVICE 
typename BackendT::UInt32_v Philox<BackendT>::Mulhilo32(typename BackendT::UInt32_v a,
							typename BackendT::UInt32_v b,
							typename BackendT::UInt32_v *hip)
{
  using UInt32_v = typename BackendT::UInt32_v;                              
  using UInt64_v = typename BackendT::UInt64_v;                              

  UInt64_v product;
  UInt32_v result;

  // UInt64_v product = ((UInt64_v)a)*((UInt64_v)b);        
  // *hip = product >> 32;                                  
  // return (UInt32_v)product; 

  // no conversion between UInt64_v and UInt32_v: size(UInt32_v)/size(UInt64_v) = 2
  // use the scalar-loop: vector->scalar bottleneck 

  size_t ibase = 0;
  size_t vsize64 = VectorSize<UInt64_v>();

  for(size_t i = 0 ; i < VectorSize<UInt32_v>()/vsize64 ; ++i) { 
    for(size_t j = 0 ; j < vsize64 ; ++j) { 
      product[j] = ((uint64_t)a[ibase+j])*((uint64_t)b[ibase+j]);
      (*hip)[ibase+j] = product[j] >> 32;                                  
      result[ibase+j] = ((uint32_t)product[j]);
    }
    ibase += vsize64;
  }
  return result; 
}

template <>
inline
VECCORE_ATT_HOST_DEVICE 
ScalarBackend::UInt32_v Philox<ScalarBackend>::Mulhilo32(uint32_t a, uint32_t b, uint32_t *hip)
{
#ifndef __CUDA_ARCH__
  uint64_t product = ((uint64_t)a)*((uint64_t)b);
  *hip = product >> 32;                                  
  return (uint32_t)product; 
#else
  *hip = __umulhi(a,b);
  return a*b;
#endif
}

template <class BackendT>
inline
VECCORE_ATT_HOST_DEVICE void Philox<BackendT>::Philox4x32bumpkey(R123::array_t<BackendT,2> key)
{
  using UInt32_v = typename BackendT::UInt32_v;                              
  key[0] += (UInt32_v)R123::PHILOX_N4_W32_0;                                       
  key[1] += (UInt32_v)R123::PHILOX_N4_W32_1;                                       
}  

template <class BackendT>
inline
VECCORE_ATT_HOST_DEVICE void 
Philox<BackendT>::Philox4x32round(R123::array_t<BackendT,4> ctr, 
                                  R123::array_t<BackendT,2> key)
{
  using UInt32_v = typename BackendT::UInt32_v;                              

  UInt32_v hi0;                                                        
  UInt32_v hi1;                                                        

  UInt32_v lo0 = Mulhilo32((UInt32_v)R123::PHILOX_M4_W32_0, ctr[0], &hi0);            
  UInt32_v lo1 = Mulhilo32((UInt32_v)R123::PHILOX_M4_W32_1, ctr[2], &hi1);            

  // a significant performance bottleneck for evaluating ctr[0] for both 
  // scalar and vector backends, especially in the assignment op= (to hi1)
  ctr[0] = hi1^ctr[1]^key[0]; 
  ctr[1] = lo1;
  ctr[2] = hi0^ctr[3]^key[1];
  ctr[3] = lo0;
}  

template <class BackendT>
VECCORE_ATT_HOST_DEVICE 
void Philox<BackendT>::Gen(R123::array_t<BackendT,4> ctr,
                           R123::array_t<BackendT,2> key,
                           R123::array_t<BackendT,4> output)
{
  // 10 rounds
  Philox4x32round(ctr, key); Philox4x32bumpkey(key);   
  Philox4x32round(ctr, key); Philox4x32bumpkey(key);  
  Philox4x32round(ctr, key); Philox4x32bumpkey(key);  
  Philox4x32round(ctr, key); Philox4x32bumpkey(key);  
  Philox4x32round(ctr, key); Philox4x32bumpkey(key);  
  Philox4x32round(ctr, key); Philox4x32bumpkey(key);  
  Philox4x32round(ctr, key); Philox4x32bumpkey(key);  
  Philox4x32round(ctr, key); Philox4x32bumpkey(key);  
  Philox4x32round(ctr, key); Philox4x32bumpkey(key);  
  Philox4x32round(ctr, key);

  //  return ctr;                                                         
  for (int i = 0 ; i < 4 ; i++) { 
    output[i] = ctr[i];                                                         
  }
}

//Copy a scalar state to at the ith-lane of the vector state
template <>
VECCORE_FORCE_INLINE
VECCORE_ATT_HOST_DEVICE
void Philox<VectorBackend>::SetStateAt(unsigned int i, State_s *state) {
  for(int j=0 ; j < 2 ; ++j) {
    this->fState->key[j][i] = state->key[j];
  }
  for(int j=0 ; j < 4 ; ++j) {
    this->fState->ctr[j][i] = state->key[j];
    this->fState->ukey[j][i] = state->key[j];
  }
}

//Return the ith-lane of the vector state to a scalar state
template <>
VECCORE_FORCE_INLINE
VECCORE_ATT_HOST_DEVICE
Philox<VectorBackend>::State_s Philox<VectorBackend>::GetStateAt(unsigned int i) 
{ 
  State_s state;

  for(int j = 0 ; j < 2 ; ++j) {
    state.key[j]= this->fState->key[j][i];
  }
  for(int j = 0 ; j < 4 ; ++j) {
    state.ctr[j]= this->fState->ctr[j][i];
    state.ukey[j]= this->fState->ukey[j][i];
  }

  return state;
} 
 
} // end namespace impl
} // end namespace vecRng

#endif
