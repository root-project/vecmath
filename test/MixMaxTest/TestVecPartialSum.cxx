// Test Partial Sum

#include <iostream>
#include <iomanip>
#include "timer.h"

// using dataType = long int; 

// #include "PartialSumVec.h"
#include "VecMath/Rng/PartialSumVec.h"

constexpr unsigned int Nscalar=16, Nvec= (Nscalar+3)/4;

using longType = UInt64_v;   // unsigned long int;
using vecType = VcSIMDArray<4>::UInt64_v;

template <typename uLongType, int num>
void SimplePartialSum(const uLongType inArr[], uLongType outArr[], uLongType overflowArr[] = 0 )
{
  outArr[0] = inArr[0];
  uLongType  overflow= 0;
  for( int i= 1; i< num; i++){
    outArr[i] = outArr[i-1] + inArr[i];
    if( outArr[i] < outArr[i-1] )
       overflow++;
    if( overflowArr )
       overflowArr[i] = overflow;
  }
}

int TestPartialSum()
{
  longType arrA[Nscalar], arrPsum[Nscalar], overflowArr[Nscalar], outPSnew[Nscalar];
  vecType  vecInpt[Nvec], vecOutp[Nvec];
  
  for( int i= 0; i< Nscalar; i++) {
    arrA[i]    = i;
    arrPsum[i] = 0;
  }
  for( int ind= 0; ind*4 < Nvec; ind++) {
    int i= 4 * ind;
    vecInpt[i]    = { arrA[i],
                      i+1<Nscalar ? arrA[i+1] : 0UL,
                      i+2<Nscalar ? arrA[i+2] : 0UL,
                      i+3<Nscalar ? arrA[i+3] : 0UL };
    Vecoutp[i]    = { 0L, 0L, 0L,  0L };
  }
  
  SimplePartialSum<longType, Nscalar> ( arrA, /* Nscalar,*/ arrPsum, overflowArr );

  // Now do an in-place partial sum
  PartialSumV1<vecType, Nvec>( vecInpt, vecOutp );

  // New partial sum
  // PartialSumV2<longType, Nscalar>( arrA, outPSnew );  

  // Check
  int nBad= 0;
  for( int i= 0; i< Nscalar; i++)
  {
    longType vecRes= vecCore::Get( vecOutp[i/4], i%4 );
    if( arrA[i] != vecRes ) {
      nBad ++; 
      std::cerr << " Error in location " << i 
               << "  scalar value  = " << arrA[i] 
               << "  vector result = " << vecRes << std::endl;
      // fprintf( stderr, " Error in location %d : value= %ld  - expected %ld \n",
      //               i, arrA[i], arrPsum[i] );
    }
  }

  return nBad;
}

/********
int BenchmarkPartialSums(const unsigned int nRepetitions= 1000 )
{
  longType arrA[Nscalar], arrOutNew[Nscalar], arrOutOld[Nscalar];
  static Timer<nanoseconds> timer;

  // longType result[nRepetitions];
  // longType resold[nRepetitions];

  longType* result= new longType[nRepetitions];
  longType* resold= new longType[nRepetitions];
  
  for( int i= 0; i< nRepetitions; i++)
    result[i]= 0;

  for( int i= 0; i< Nscalar; i++){
    arrA[i]= i;
    arrOutOld[i]= 0;
  }
  
  timer.Start();
  for( int j= 0; j< nRepetitions; j++)  
  {
     for( int i= 0; i< Nscalar; i++)
        arrA[i]= i;
     arrA[ j % Nscalar] = 0;

     PartialSumV1<longType, Nscalar>( arrA );
     result[j]= arrA[Nscalar-1];
  }
  double newtime = timer.Elapsed();

  // Store result of last partial sums
  for( int i= 0; i< Nscalar; i++)
     arrOutNew[i] = arrA[i];
     
  timer.Start();
  for( int j= 0; j < nRepetitions; j++)  
  {
     for( int i= 0; i< Nscalar; i++)
        arrA[i]= i;
     arrA[ j % Nscalar] = 0;

     SimplePartialSum( arrA, Nscalar, arrOutOld );
     resold[j]= arrOutOld[Nscalar-1];
  }
  double oldtime = timer.Elapsed();  


  // std::cout << " Results for N= " << Nscalar << std::endl;
  printf( " Results for N= %4d  - number of repetitions=  %8d K \n", Nscalar, nRepetitions / 1000);
  
  // Check last result
  int nBad= 0;
  for( int i= 0; i< Nscalar; i++)
  {
    if( arrOutNew[i] != arrOutOld[i] ) { 
      nBad ++; 
      std::cerr << " Last result> Error in location " << i 
                << "  value = " << arrOutNew[i]
                << "  expected = " << arrOutOld[i] << std::endl;
      // fprintf( stderr, " Error in location %d : value= %ld  - expected %d \n",
      //                 i, arrA[i], i+1 );
    }
  }

  // Check all sums
  int nBad2= 0;
  for( int i= 0; i< nRepetitions; i++)
  {
    if( result[i] != resold[i] ) { 
      nBad2 ++; 
      // std::cerr << " Error in location " << i 
      //       << "  value = " << arrA[i] 
      //       << "  expected = " << i+1 << std::endl;
      fprintf( stderr, " Error in sum of repetition  %d : value= %ld  - expected %ld \n",
                       i, result[i], resold[i] );
    }
  }

//  std::cout << "Time of templated method= " << std::setprecision(3) << newtime / nRepetitions << " ns "
//            << " (per trial ) " << std::endl
//            << "     of simple    method= " << oldtime / nRepetitions << " ns " << std::endl;

  fprintf( stdout, "Time of templated method= %5.1f ns (per trial)\n", newtime / nRepetitions);
  fprintf( stdout, "     of   simple  method= %5.1f ns (per trial)\n", oldtime / nRepetitions);

  delete[] result;
  delete[] resold; 
  
  return nBad + nBad2;
}
 ***/
int main() // int argc, char **argv)
{
  // if( argc == 1) N = atoi( argv[0] );
  TestPartialSum();

  //  int nRepetitions= 1000000;
  //  BenchmarkPartialSums(nRepetitions);
}
