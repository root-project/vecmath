#include "RngBenchmarker.h"
#include "RngBenchmarker_cpu.h"

#include "VecCore/Timer.h"

namespace vecrng {

RngBenchmarker::RngBenchmarker()
  : fNSample(100), fRepetition(1)
{
}

RngBenchmarker::~RngBenchmarker()
{
}

int RngBenchmarker::RunBenchmark()
{
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

int RngBenchmarker::RunBenchmarkRng()
{
  if (fVerbosity > 0) {
  }

  RunTest();
  std::cout << " ---- Running Scalar()\n";  
  RunScalar();
  std::cout << " ---- Running Vector()  - ie. type 1 \n";
  RunVector();
  std::cout << " ---- Running Vector2() - ie. type 2 \n";  
  RunVector2();
  std::cout << " ---- Running NState()\n";  
  RunNState();

#ifdef RNGTEST_MKL
  RunMKLVSL();
#endif

#ifdef RNGTEST_CUDA
    RunCuda();
#endif

  return 0;
}

double variance( double sumValues, double sumSquares, int repetitions )
{
   double var=0.0;
   double meanTime= sumValues / repetitions;
   double varSq= sumSquares / repetitions - meanTime * meanTime;
   if( varSq > 0.0 ) var = std::sqrt( varSq );
   return var;
}

void PrintReport(const char* method, const char* version, double meanTime,
                 double sigmaTime, double resultTotal )
{
  printf(" %-15s  %-13s Time = ", method, version);
      
  if( meanTime*1E-6 > 1.0 ) {
     const double msMult = 1.e-6; // micro-sec per ns      
     printf("%6.3f +- %6.3f msec ", meanTime*msMult, sigmaTime*msMult);
  } else if( meanTime*1E-3 > 1.0 ) {     
     const double usMult = 1.e-3; // micro-sec per ns 
     printf("%6.3f +- %6.3f usec ", meanTime*usMult, sigmaTime*usMult );
  } else {
     const double nsMult = 1.0; //  directly ns
     printf("%6.1f +- %6.1f nsec ", meanTime*nsMult, sigmaTime*nsMult );
  }
  printf( "Sum = %g\n", resultTotal);  
}

void ProcessAndPrint( const char* method, const char* version,
                      double elapsedTotal,
                      double elapsedTotSq,
                      int    nRepetition,
                      int    nSamples,
                      double resultTotal )
{
  double meanTime = elapsedTotal/nRepetition;
  double sigmaTime= variance( elapsedTotal, elapsedTotSq, nRepetition );
  // meanTime  /= nSamples;
  // sigmaTime /= nSamples;  
  PrintReport(method, version, meanTime, sigmaTime, resultTotal );
}

// Reference
void RngBenchmarker::RunTest()
{
  // test GNU rand() as a reference

  double elapsedTotal = 0.;
  double elapsedTotSq = 0.;
  double resultTotal = 0.;
#ifdef OLDVAR
  double *trialTime = new double [fRepetition];;
#endif
  static Timer<nanoseconds> timer;

  for (unsigned r = 0; r < fRepetition; ++r) {

    double result = 0.;

    timer.Start();

    for (int i = 0; i < fNSample ; ++i) {
      result += (double)rand()/RAND_MAX;
    }
    double time = timer.Elapsed();
    elapsedTotal += time;
    elapsedTotSq += time * time;
    resultTotal += result;
  }

#if 1  
  ProcessAndPrint( "TestRand", "std::rand()", elapsedTotal, elapsedTotSq, fRepetition, fNSample,
                   resultTotal );  
#else
  double meanTime = elapsedTotal/fRepetition;
  double sigmaTime= variance( elapsedTotal, elapsedTotSq, fRepetition );
  // PrintReport( "TestRand", "std::rand()", meanTime, sigmaTime, resultTotal );
  
  printf(" %-15s std::rand()   ", "TestRand");
  if( meanTime*1E-6 > 0.1 ) {
    printf("Time = %4.3f +- %4.3f msec ", // \n", 
  	    meanTime*1E-6, sigmaTime*1E-6);
  } else {
    const double usMult = 1.e-3; // micro-sec per ns 
    printf( "Time= %4.3f +- %4.3f usec ",
            meanTime*usMult, sigmaTime*usMult );     
  }
  printf( "Sum = %g\n", resultTotal);  
//  printf(" Check value= %4.3f", sqrt(sumDiff/fRepetition) );
  printf("\n");
#endif  
}

// Scalar
void RngBenchmarker::RunScalar()
{
  double meanTime[kNumberRng];
  double sigmaTime[kNumberRng];
  double resultTotal[kNumberRng];
  double varncTime[kNumberRng];
  
  double result = 0;
  // double *trialTime = new double [fRepetition];;

  for (unsigned int k = 0; k < kNumberRng ; ++k) {

    meanTime[k] = 0.;
    sigmaTime[k] = 0.;
    varncTime[k] = 0.;
    resultTotal[k] = 0.;

    double elapsedTotal = 0.0, elapsedTotalSq = 0.0;

    for (unsigned r = 0; r < fRepetition; ++r) {
       // trialTime[r] = 0.0;
      result = 0.0;
      double time = ScalarKernelFunc[k](fNSample,result);
      elapsedTotal += time;
      elapsedTotalSq += time * time;
      resultTotal[k] += result;
    }

    // meanTime[k] = elapsedTotal/fRepetition;
    // sigmaTime[k] = variance( elapsedTotal, elapsedTotalSq, fRepetition );       
    // varncTime[k] = sqrt( elapsedTotalSq/fRepetition - meanTime[k] * meanTime[k] );

    ProcessAndPrint( RngName[k], "ScalarBackend", elapsedTotal, elapsedTotalSq,
                     fRepetition, fNSample, resultTotal[k] );
  }
}

// Vector 
void RngBenchmarker::RunVector()
{
  double meanTime[kNumberVectorRng];
  double sigmaTime[kNumberVectorRng];
  double resultTotal[kNumberVectorRng];

  double result = 0;
  // double *trialTime = new double [fRepetition];;

  for (unsigned int k = 0; k < kNumberVectorRng ; ++k) {

    meanTime[k] = 0.;
    sigmaTime[k] = 0.;
    resultTotal[k] = 0.;

    //std::cout << "RunVector> k= " << k << " VecRNG is " << VecRngName[k] << std::endl;
    
    double elapsedTotal = 0.0, elapsedTotalSq = 0.0;
    
    for (unsigned r = 0; r < fRepetition; ++r) {
      // trialTime[r] = 0.0;
      result = 0.0;
      double time = VectorKernelFunc[k](fNSample,result);
      elapsedTotal += time;
      elapsedTotalSq += time * time;
      resultTotal[k] += result;
    }
#if 1
    ProcessAndPrint( VecRngName[k], "VectorBackend", elapsedTotal, elapsedTotalSq,
                     fRepetition, fNSample, resultTotal[k] );
#else    
    meanTime[k] = elapsedTotal/fRepetition;
    sigmaTime[k] = variance( elapsedTotal, elapsedTotalSq, fRepetition );
    
    printf(" %-15s  VectorBackend Time = ", VecRngName[k] );
    if( meanTime[k]*1E-6 > 0.1 ) {
      const double msMult = 1.e-6; // millisecond per nanosec
      printf( "%4.3f +- %4.3f msec ", 
              meanTime[k]*msMult, sigmaTime[k]*msMult );
    } else {
      const double usMult = 1.e-3; // micro-sec per ns 
      printf( "%4.3f +- %4.3f usec ",
              meanTime[k]*usMult, sigmaTime[k]*usMult );
    }
    printf( "Sum = %g\n", resultTotal[k]);
#endif    
  }
}

// Vector 
void RngBenchmarker::RunVector2()
{
  double meanTime[numVecKnlFuncs];
  double sigmaTime[numVecKnlFuncs];
  double resultTotal[numVecKnlFuncs];
  double result = 0;

  for (unsigned int k = 0; k < numVecKnlFuncs ; ++k) {

    meanTime[k] = 0.;
    sigmaTime[k] = 0.;
    resultTotal[k] = 0.;

    double elapsedTotal = 0.0, elapsedTotSq = 0.0;
    auto funcAndName = VectorKnlFuncAndName[k];
    auto func= std::get<0>(funcAndName); // first();
    // auto name=  std::get<std::string>(funcAndName);  // funcAndName.second();
    for (unsigned r = 0; r < fRepetition; ++r) {
      result = 0.0;
      // trialTime[r] = VectorKernelFunc[k](fNSample,result);
      // trialTime[r] = func(fNSample,result);
      double oneTm = func(fNSample,result);
      elapsedTotal += oneTm; // trialTime[r];
      elapsedTotSq += oneTm * oneTm;
      resultTotal[k] += result;
    }

    std::string funcName=  std::get<1 /*std::string*/ >(funcAndName);  // funcAndName.second();    
    ProcessAndPrint( funcName.c_str(), "VectorBackend", elapsedTotal, elapsedTotSq,
                     fRepetition, fNSample, resultTotal[k] );
  }
}

   
void RngBenchmarker::RunNState()
{
  double meanTime[kNumberRng];
  double sigmaTime[kNumberRng];
  double resultTotal[kNumberRng];

  double result = 0;
  // double *trialTime = new double [fRepetition];
  double trialTime;

  for (unsigned int k = 0; k < kNumberRng ; ++k) {

    meanTime[k] = 0.;
    sigmaTime[k] = 0.;
    resultTotal[k] = 0.;

    double elapsedTotal = 0.0, elapsedTotSq = 0.0;

    for (unsigned r = 0; r < fRepetition; ++r) {
       // trialTime[r] = 0.0;
      result = 0.0;
      trialTime= StateKernelFunc[k](fNSample,result);
      elapsedTotal += trialTime;
      elapsedTotSq += trialTime * trialTime;
      resultTotal[k] += result;
    }
#if 1    
    ProcessAndPrint( RngName[k], "ScalarNstates", elapsedTotal, elapsedTotSq,
                     fRepetition, fNSample, resultTotal[k] );
#else    
    meanTime[k] = elapsedTotal/fRepetition;
    sigmaTime[k]= variance( elapsedTotal, elapsedTotSq, fRepetition );
  }
  // delete [] trialTime;

  for (int k = 0; k < kNumberRng; ++k) {
    if( meanTime[k]*1E-6 > 0.1 ) {     
       printf(" %-15s  ScalarNstates Time = %4.3f +- %4.3f msec Sum = %g\n", 
              RngName[k], meanTime[k]*1E-6, sigmaTime[k]*1E-6, resultTotal[k]);
    } else {
       printf(" %-15s  ScalarNstates Time = %4.3f +- %4.3f usec Sum = %g\n", 
              RngName[k], meanTime[k]*1E-3, sigmaTime[k]*1E-3, resultTotal[k]);
    }
#endif      
  }
}

#ifdef RNGTEST_MKL
void RngBenchmarker::RunMKLVSL()
{
  double meanTime[kNumberVsl];
  double sigmaTime[kNumberVsl];
  double resultTotal[kNumberVsl];

  double result =0;
  double *trialTime = new double [fRepetition];;

  for (unsigned int k = 0; k < kNumberVsl ; ++k) {

    meanTime[k] = 0.;
    sigmaTime[k] = 0.;
    resultTotal[k] = 0.;

    double elapsedTotal = 0.0;

    for (unsigned r = 0; r < fRepetition; ++r) {
      trialTime[r] = 0.0;
      result = 0.0;
      trialTime[r] = VSLKernelFunc[k](fNSample,result);
      elapsedTotal += trialTime[r];
      resultTotal[k] += result;
    }

    meanTime[k] = elapsedTotal/fRepetition;
    double sumDiff = 0;

    for (unsigned r = 0; r < fRepetition; ++r) {
      double delta  = (trialTime[r] - meanTime[k]);
      sumDiff += delta*delta;
    }
    sigmaTime[k] = sqrt(sumDiff/fRepetition);   
  }
  delete [] trialTime;

  for (int k = 0; k < kNumberVsl; ++k) {
    printf(" %-15s  Intel-MKL/VSL Time = %3.3f +- %3.3f msec Sum = %g\n", 
    VslName[k], meanTime[k]*1E-6, sigmaTime[k]*1E-6, resultTotal[k]);
  }
}
#endif

} // end namespace vecrng
