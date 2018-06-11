#include "VecCore/Timer.h"

#include "VecMath/Rng/MRG32k3a.h"
#include "VecMath/Rng/Threefry.h"
#include "VecMath/Rng/Philox.h"

#include "VecCoreLib/Rng/JoiningProxyVecRNG.h"

// #define NEW_TEST_JOINING_PROXY    1

#ifdef NEW_TEST_JOINING_PROXY
#include "VecCoreLib/Rng/JoiningProxyVecRNG.h"
#endif

#ifdef RNGTEST_MKL
#include "mkl_vsl.h"
#include <omp.h>
#endif

#include "RngTest.h"

namespace vecrng {

// Scalar

double ScalarMRG32k3a(int nsample, double& result)
{
  // Scalar MRG32k3a
  static vecRng::cxx::MRG32k3a<ScalarBackend> rng;
  rng.Initialize();

  static Timer<nanoseconds> timer;
  double elapsedTime = 0.;

  double sum = 0;

  timer.Start();

  for (int i = 0; i < nsample ; ++i) {
    sum += rng.Uniform<ScalarBackend>();
  }

  elapsedTime = timer.Elapsed();
  result = sum;

  return elapsedTime;
}
   
// ---------------============================------------------------------------------

double ScalarThreefry(int nsample, double& result)
{
  // Scalar Threefry
  static vecRng::cxx::Threefry<ScalarBackend> rng;
  rng.Initialize();

  static Timer<nanoseconds> timer;
  double elapsedTime = 0.;

  double sum = 0;

  timer.Start();

  for (int i = 0; i < nsample ; ++i) {
    sum += rng.Uniform<ScalarBackend>();
  }

  elapsedTime = timer.Elapsed();
  result = sum;

  return elapsedTime;
}

// ---------------============================------------------------------------------

double ScalarPhilox(int nsample, double& result)
{
  // Scalar Philox - use scalar method of VectorBackend (temporarily)

  static vecRng::cxx::Philox<ScalarBackend> rng;
  rng.Initialize();

  static Timer<nanoseconds> timer;
  double elapsedTime = 0.;

  double sum = 0;

  timer.Start();

  for (int i = 0; i < nsample ; ++i) {
    sum += rng.Uniform<ScalarBackend>();
  }

  elapsedTime = timer.Elapsed();
  result = sum;

  vecRng::cxx::Philox<ScalarBackend> rng1 = rng;

  return elapsedTime;
}

// ---------------============================------------------------------------------
// Vector
// ---------------============================------------------------------------------

// #ifndef ORIGINAL_VECTOR_MRG
double VectorMRG32k3a_JoinByHand(int nsample, double& result)
{
  // Vector MRG32k3a
  using Double_v = typename VectorBackend::Double_v;
  constexpr int vsize = VectorSize<Double_v>();

  // Multiple scalar MRG32k3a
  using ScalarRngType = vecRng::cxx::MRG32k3a<ScalarBackend>;

  // Multiple scalar MRG32k3a
  static ScalarRngType // vecRng::cxx::MRG32k3a<ScalarBackend>
          scalarRng[vsize];
  for (int i = 0; i < vsize ; ++i) {  
     scalarRng[i].Initialize(i);
  }

  vecRng::cxx::MRG32k3a<VectorBackend> rng;
  rng.Initialize();

  // Stuff the (vector) rng with the contents of the scalar ones !!
   for( int i= 0; i < vsize; ++i )
      rng.SetPartOfState( scalarRng[i].GetState(), i);
  
  static Timer<nanoseconds> timer;
  double elapsedTime = 0.;

  Double_v sum = 0.;
  int      ntrials=  nsample/vsize;
  timer.Start();

  for (int i = 0; i < ntrials ; ++i) {
    sum += rng.Uniform<VectorBackend>();
  }

  elapsedTime = timer.Elapsed();
  for (int i = 0; i < vsize ; ++i) result += sum[i];

  using std::cout; using std::endl;
  cout << "Next values: " << endl;
  for (int i = 0; i < vsize ; ++i) 
     cout << " [ " << i << " ] " << scalarRng[i].Uniform<ScalarBackend>() << endl;
  
  return elapsedTime;
}
// #else
// ---------------============================------------------------------------------
double VectorMRG32k3a_Basic(int nsample, double& result)
{
  // Vector MRG32k3a
  using Double_v = typename VectorBackend::Double_v;
  int vsize = VectorSize<Double_v>();

  vecRng::cxx::MRG32k3a<VectorBackend> rng;
  rng.Initialize();

  static Timer<nanoseconds> timer;
  double elapsedTime = 0.;

  Double_v sum = 0.;

  timer.Start();

  for (int i = 0; i < nsample/vsize ; ++i) {
    sum += rng.Uniform<VectorBackend>();
  }

  elapsedTime = timer.Elapsed();
  for (int i = 0; i < vsize ; ++i) result += sum[i];

  return elapsedTime;
}
// #endif

double VectorMRG32k3a_JoiningProxy(int nsample, double& result)
{
   // Vector MRG32k3a
  using Double_v = typename VectorBackend::Double_v;
  constexpr int vsize = VectorSize<Double_v>();

  using ScalarRngType = vecRng::cxx::MRG32k3a<ScalarBackend>;
  
  static ScalarRngType scalarRng[vsize];
  static ScalarRngType * arrayScalarRngPtr[vsize];
  for (int i = 0; i < vsize ; ++i) {  
     scalarRng[i].Initialize(i);
     arrayScalarRngPtr= & scalarRng[i];
  }
  using VectorRngProxyType = JoiningProxyVecRNG< vecRng::cxx::MRG32k3a<VectorBackend>;

  // Create a scope, so see also the copying back of scalar RNGs
  { 
     // Prepare the (vector) proxy rng giving the scalar ones !!
     VectorRngProxyType rng( arrayScalarRngPtr, vsize );
     
     static Timer<nanoseconds> timer;
     double elapsedTime = 0.;
     
     Double_v sum = 0.;
     int      ntrials=  nsample/vsize;
     timer.Start();
     
     for (int i = 0; i < ntrials ; ++i) {
        sum += rng.Uniform<VectorBackend>();
     }

     elapsedTime = timer.Elapsed();
     for (int i = 0; i < vsize ; ++i) result += sum[i];
  }

  using std::cout;
  using std::endl;  
  cout << "Next values: " << endl;
  for (int i = 0; i < vsize ; ++i) 
     cout << " [ " << i << " ] " << scalarRng[i].uniform() << endl;
     
  return elapsedTime;
}
double VectorMRG32k3a(int nsample, double& result)
{
   double time=0.0;
   time=
      // VectorMRG32k3a_Basic(nsample, result);
      // VectorMRG32k3a_JoinByHand(nsample, result);      
         VectorMRG32k3a_JoiningProxy(nsample, result);   
   
   return time;
}
// ---------------============================------------------------------------------

#ifdef NEW_TEST_JOINING_PROXY
double VectorJoiningMRG32k3a(int nsample, double& result)
{
  // Vector MRG32k3a
  using Double_v = typename VectorBackend::Double_v;
  // using template MRG32k3a<Backend> = typename vecRng::cxx::MRG32k3a<Backend>;  
  using JoiningProxyMRG32K3a = typename vecRng::cxx::JoiningProxyVecRNG<MRG32k3a<VectorBackend>,
                                                                        MRG32k3a<ScalarBackend>,
                                                                        VectorBackend>;
  constexpr int VecSize = VectorSize<Double_v>();

  // Multiple scalar MRG32k3a
  static vecRng::cxx::MRG32k3a<ScalarBackend> scalarRng[VecSize];
  for (int i = 0; i < VecSize ; ++i) {  
     scalarRng[i].Initialize(i);
  }
  
  // vecRng::cxx::MRG32k3a<VectorBackend> rng;
  // rng.Initialize();

  JoiningProxyMRG32k3a rng( scalarRng );
  // rng.Join( scalarRng );
  
  static Timer<nanoseconds> timer;
  double elapsedTime = 0.;

  Double_v sum = 0.;

  timer.Start();

  for (int i = 0; i < nsample/VecSize ; ++i) {
    sum += rng.Uniform<VectorBackend>();
  }

  elapsedTime = timer.Elapsed();
  for (int i = 0; i < VecSize ; ++i) result += sum[i];

  // rng.Split();
  
  return elapsedTime;
}
#endif

// ---------------============================------------------------------------------
   
double VectorThreefry(int nsample, double& result)
{
  // Vector Threefry
  using Double_v = typename VectorBackend::Double_v;
  int vsize = VectorSize<Double_v>();

  static vecRng::cxx::Threefry<VectorBackend> rng;
  rng.Initialize();

  static Timer<nanoseconds> timer;
  double elapsedTime = 0.;

  Double_v sum = 0;

  timer.Start();

  for (int i = 0; i < nsample/vsize ; ++i) {
    sum += rng.Uniform<VectorBackend>();
  }

  elapsedTime = timer.Elapsed();
  for (int i = 0; i < vsize ; ++i) result += sum[i];

  return elapsedTime;
}

// ---------------============================------------------------------------------

double VectorPhilox(int nsample, double& result)
{
  // Vector Philox
  using Double_v = typename VectorBackend::Double_v;
  int vsize = VectorSize<Double_v>();

  vecRng::cxx::Philox<VectorBackend> rng;
  rng.Initialize();

  static Timer<nanoseconds> timer;
  double elapsedTime = 0.;

  Double_v sum = 0;

  timer.Start();

  for (int i = 0; i < nsample/vsize ; ++i) {
    sum += rng.Uniform<VectorBackend>();
  }

  elapsedTime = timer.Elapsed();
  for (int i = 0; i < vsize ; ++i) result += sum[i];

  return elapsedTime;
}

// ---------------============================------------------------------------------

double StateMRG32k3a(int nsample, double& result)
{
  // Scalar MRG32k3a
  static vecRng::cxx::MRG32k3a<ScalarBackend> rng;

  int theNBlocks = 26;
  int theNThreads = 192;

  vecRng::MRG32k3a<ScalarBackend>::State_t* hstates 
    = (vecRng::MRG32k3a<ScalarBackend>::State_t *) malloc (theNBlocks*theNThreads*sizeof(vecRng::MRG32k3a<ScalarBackend>::State_t));
  rng.Initialize(hstates,theNBlocks*theNThreads);

  static Timer<nanoseconds> timer;
  double elapsedTime = 0.;

  double sum = 0;

  timer.Start();

  for(int i = 0 ; i < theNBlocks*theNThreads ; ++i) {
    for (int j = 0; j < nsample/(theNBlocks*theNThreads) ; ++j) {
      sum += rng.Uniform<ScalarBackend>(&hstates[i]);
    }
  }

  elapsedTime = timer.Elapsed();
  result = sum;

  return elapsedTime;
}

// ---------------============================------------------------------------------

double StateThreefry(int nsample, double& result) 
{
  // Scalar Threefry
  static vecRng::cxx::Threefry<ScalarBackend> rng;

  int theNBlocks = 26;
  int theNThreads = 192;

  vecRng::Threefry<ScalarBackend>::State_t* hstates 
    = (vecRng::Threefry<ScalarBackend>::State_t *) malloc (theNBlocks*theNThreads*sizeof(vecRng::Threefry<ScalarBackend>::State_t));
  rng.Initialize(hstates,theNBlocks*theNThreads);

  static Timer<nanoseconds> timer;
  double elapsedTime = 0.;

  double sum = 0;

  timer.Start();

  for(int i = 0 ; i < theNBlocks*theNThreads ; ++i) {
    for (int j = 0; j < nsample/(theNBlocks*theNThreads) ; ++j) {
      sum += rng.Uniform<ScalarBackend>(&hstates[i]);
    }
  }

  elapsedTime = timer.Elapsed();
  result = sum;

  return elapsedTime;
}

// ---------------============================------------------------------------------

double StatePhilox(int nsample, double& result) 
{
  // Scalar Philox
  //  static vecrng::cxx::VecRNG<Philox<ScalarBackend>, ScalarBackend, Philox_t<ScalarBackend>> rng;
  static vecRng::cxx::Philox<ScalarBackend> rng;

  int theNBlocks = 26;
  int theNThreads = 192;

  vecRng::Philox<ScalarBackend>::State_t* hstates
    = (vecRng::Philox<ScalarBackend>::State_t *) malloc (theNBlocks*theNThreads*sizeof(vecRng::Philox<ScalarBackend>::State_t));
  rng.Initialize(hstates,theNBlocks*theNThreads);

  static Timer<nanoseconds> timer;
  double elapsedTime = 0.;

  double sum = 0;

  timer.Start();

  for(int i = 0 ; i < theNBlocks*theNThreads ; ++i) {
    for (int j = 0; j < nsample/(theNBlocks*theNThreads) ; ++j) {
      sum += rng.Uniform<ScalarBackend>(&hstates[i]);
    }
  }

  elapsedTime = timer.Elapsed();
  result = sum;

  return elapsedTime;
}

// ---------------============================------------------------------------------

// Intel VMK/VSL
#ifdef RNGTEST_MKL

double VSLRngTest(int nsample, double& result, int rngType)
{
  static Timer<nanoseconds> timer;
  double elapsedTime = 0.;

  //buffer size
  const int N = 32;

  double r[N];
  double a=0.0;
  double b=1.0;

  double sum = 0;

  // Initialize
  unsigned int seed = 7777777;
  VSLStreamStatePtr stream;
  vslNewStream(&stream, rngType, (MKL_INT)seed);

  timer.Start();

  // Call RNG
  for (int i = 0; i < nsample/N ; ++i) {
    vdRngUniform( VSL_RNG_METHOD_UNIFORM_STD_ACCURATE, stream, N, r, a, b);
    for(int j = 0 ; j < N ; ++j) sum += r[j];
  }
  elapsedTime = timer.Elapsed();

  result = sum;

  // Deinitialize
  vslDeleteStream( &stream );

  return elapsedTime;
}

#endif

} // end namespace vecrng
