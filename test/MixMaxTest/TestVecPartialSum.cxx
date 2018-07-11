// Test Partial Sum

#include <iostream>
#include <iomanip>
#include <cassert>

#include "timer.h"

// using dataType = long int; 

#include "VecCore/VecCore"
#include "VecCore/backend/UMESimdArray.h"

#include "VecMath/Rng/PartialSumVec.h"

constexpr unsigned int Nscalar=16, Nvec= (Nscalar+3)/4;

using scalarLong = vecCore::UInt64_s;   // unsigned long int;
// using vecType = VcSIMDArray<4>::UInt64_v;
using vecType = vecCore::backend::UMESimdArray<4>::UInt64_v;

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
  scalarLong arrA[Nscalar], arrPsum[Nscalar], overflowArr[Nscalar], outPSnew[Nscalar];
  vecType  vecInpt[Nvec], vecOutp[Nvec], tokenVT;
  // const vecCore::Size_s vecWidth= vecCore::VectorSize<decltype(tokenVT)>;
  constexpr int vecWidth= vecCore::VectorSize<decltype(tokenVT)>();
  assert( vecWidth == 4  );

  bool verbose= true;
  if( verbose ) {
     std::cout << "Nscalar = " << Nscalar << " Nvec= " << Nvec << std::endl;
  }
  
  for( int i= 0; i< Nscalar; i++) {
    arrA[i]    = i;
    arrPsum[i] = 0;
  }

  for( int ind= 0; ind < Nvec; ind++) {
    int i= 4 * ind;
    vecInpt[ind]  = { arrA[i],
                      i+1<Nscalar ? arrA[i+1] : 0UL,
                      i+2<Nscalar ? arrA[i+2] : 0UL,
                      i+3<Nscalar ? arrA[i+3] : 0UL };
    vecOutp[ind]    = { 0L, 0L, 0L,  0L };
    
    if( verbose )  {
       std::cout << " [ " << ind << " => " << i << " ] = " ;
       for( int j= 0; j< 4; j++ ) {
          std::cout << vecCore::Get( vecInpt[ind] , j ) << " ";
       }
    }

  }
  std::cout << endl;

  SimplePartialSum<scalarLong, Nscalar> ( arrA, /* Nscalar,*/ arrPsum, overflowArr );

  // Now a vector partial sum - without modulo
  constexpr bool enableModulo = false;
  PartialSumVec<vecType, Nvec, enableModulo>( vecInpt, vecOutp);

  // New partial sum
  // PartialSumV2<longType, Nscalar>( arrA, outPSnew );  

  bool verboseOut = true;
  if( verboseOut ) { 
     // Print output - verbose mode
     printf( "   %3s    %10s   %10s  %10s \n", "I", "Input", "Scalar-psum", "Vec-par/sum");
     printf( "========================================================================\n");
     for( int i= 0; i< Nscalar; i++)
     {
        scalarLong vecRes= vecCore::Get( vecOutp[i/vecWidth], i % vecWidth );
        // std::cout << " [ " << i << " scalar value  = " << arrA[i];
        printf( " [ %3d ]  %10llu   %10llu  %10llu \n", i, arrA[i], arrPsum[i], vecRes);
     }
     printf( "========================================================================\n");
  }
  
  // Check
  int nBad= 0;
  for( int i= 0; i< Nscalar; i++)
  {
    scalarLong vecRes= vecCore::Get( vecOutp[i/vecWidth], i % vecWidth );
    if( arrPsum[i] != vecRes ) {
      nBad ++; 
      // std::cerr << " Error in location " << i 
      //          << "  scalar value  = " << arrPsum[i] 
      //          << "  vector result = " << vecRes << std::endl;
      fprintf( stderr, " Error in location %3d : scalar value= %5llu  - obtained %5llu \n",
                    i, arrPsum[i], vecRes );
    }
  }

  return nBad;
}

/********
int BenchmarkPartialSums(const unsigned int nRepetitions= 1000 )
{
  scalarLong arrA[Nscalar], arrOutNew[Nscalar], arrOutOld[Nscalar];
  static Timer<nanoseconds> timer;

  // scalarLong result[nRepetitions];
  // scalarLong resold[nRepetitions];

  scalarLong* result= new scalarLong[nRepetitions];
  scalarLong* resold= new scalarLong[nRepetitions];
  
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

     PartialSumV1<scalarLong, Nscalar>( arrA );
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
