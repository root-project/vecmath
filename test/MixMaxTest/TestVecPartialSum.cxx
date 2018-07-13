// Test Partial Sum

#include <iostream>
#include <iomanip>
#include <cassert>

#include "timer.h"

// using dataType = long int; 

// #define  __AVX2__   1   // Instruct UME::SIMD to use AVX2 ... !?
//                            Resolved by adding -march=native in CXXFLAGS
#include "VecCore/VecCore"
#include "VecCore/backend/UMESimdArray.h"

#include "VecMath/Rng/PartialSumVec.h"

//  Most important parameters
constexpr unsigned int Nscalar=16;
constexpr bool enableModulo = true;

//  Derived parameters
constexpr unsigned int Nvec= (Nscalar+3)/4;

using scalarLong = vecCore::UInt64_s;   // unsigned long int;
// using vecType = VcSIMDArray<4>::UInt64_v;
using vecType = vecCore::backend::UMESimdArray<4>::UInt64_v;
                // UME::SIMD::SIMDVec<uint64_t, 4>;
// using UInt64_v = UME::SIMD::SIMDVec<UInt64_s, N>;

// ---------------------------------------------------------------------------------

template <typename uLongType, int num>
void PartialSumWithOverflow(const uLongType inArr[], uLongType outArr[], uLongType overflowArr[] = 0 )
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

// ---------------------------------------------------------------------------------

void
PrepareInputVecArray( vecType  vecArr[Nvec],  scalarLong arrIn[Nscalar], bool verbose= false )
{
  for( int ind= 0; ind < Nvec; ind++) {
    int i= 4 * ind;
    vecArr[ind]  = { arrIn[i],
                      i+1<Nscalar ? arrIn[i+1] : 0UL,
                      i+2<Nscalar ? arrIn[i+2] : 0UL,
                      i+3<Nscalar ? arrIn[i+3] : 0UL };
    if( verbose )  {
       std::cout << " [ " << ind << " => " << i << " ] = " ;
       for( int j= 0; j< 4; j++ ) {
          std::cout << vecCore::Get( vecArr[ind] , j ) << " ";
       }
    }
  }
  if( verbose )  
    std::cout << std::endl;
    
}

// ---------------------------------------------------------------------------------

void
ConvertVecToArray( const vecType  vecArr[Nvec],
                   scalarLong     arrOut[Nscalar],
                   int  lenScalar, 
                   bool verbose= false )
{
  // auto firstElem= vecArr[0];  // Can use it to deduce vecSize
  int  vecSize = 4;   
  for( int i= 0; i < lenScalar; i++) {
     int indV = i / vecSize;
     int slot = i % vecSize;
    arrOut[i] = vecCore::Get( vecArr[indV], slot );
  }
}

void ReportSumResults( const scalarLong arrInput[Nscalar],
                       const scalarLong outScalar[Nscalar],  // Reference result
                       const vecType    vecSum[Nvec],        // New (vec) result
                       int Nscalar );

// ---------------------------------------------------------------------------------

void SimplePartialSum(const scalarLong inArr[], int num,
                      scalarLong outArr[] )
{
  outArr[0] = inArr[0];
  for( int i= 1; i< num; i++){
    outArr[i] = outArr[i-1] + inArr[i];
  }
}

// ---------------------------------------------------------------------------------

void SimplePartialSumWithModulo(const scalarLong inArr[], int num,
                                scalarLong outArr[] )
{
  outArr[0] = inArr[0];
  for( int i= 1; i< num; i++){
     outArr[i] = FastModuloMers61( outArr[i-1] + inArr[i] );
  }
}

// =================================================================================

int TestPartialSum( bool useModulo= false )
{
  scalarLong arrA[Nscalar], arrPsum[Nscalar], overflowArr[Nscalar], outPSnew[Nscalar];
  vecType  vecInpt[Nvec], vecOutp[Nvec], tokenVT;
  // const vecCore::Size_s vecWidth= vecCore::VectorSize<decltype(tokenVT)>;
  constexpr int vecWidth= vecCore::VectorSize<decltype(tokenVT)>();
  assert( vecWidth == 4  );

  bool verbose= true;
  if( verbose ) {
     std::cout << "TestPartialSum> ";
     std::cout << "Nscalar = " << Nscalar << " Nvec= " << Nvec << std::endl;
  }

  for( int i= 0; i< Nscalar; i++) {
    arrA[i]    = gMers61 - i;
    arrPsum[i] = 0;
  }

  bool verbInp= false;
  PrepareInputVecArray( vecInpt, arrA, verbInp ); 
  
  // PartialSumWithOverFlow<scalarLong, Nscalar> ( arrA, arrPsum, overflowArr );
  if( ! useModulo )
     SimplePartialSum( arrA, Nscalar, arrPsum );
  else  
     SimplePartialSumWithModulo( arrA, Nscalar, arrPsum );  

  
  // Now a vector partial sum - without modulo
  constexpr bool enabledModulo = true;

  if( ! useModulo )  
     PartialSumVec<vecType, Nvec, false >( vecInpt, vecOutp);
  else
     PartialSumVec<vecType, Nvec, enabledModulo>( vecInpt, vecOutp);     

  // New partial sum
  // PartialSumV2<longType, Nscalar>( arrA, outPSnew );  

  bool verboseOut = false;
  // Print results - by default
  
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
  
  if( verboseOut || (nBad != 0) ) {
     ReportSumResults( arrA,    // input
                       arrPsum, // reference: scalar result
                       vecOutp, // new result (vector)
                       Nscalar );
  }
  
  return nBad;
}

#define ENABLE_BENCHMARK  1


#ifdef  ENABLE_BENCHMARK
int BenchmarkPartialSums(const unsigned int nRepetitions= 1000 )
{
  scalarLong arrA[Nscalar], arrOutScalar[Nscalar], arrOutVector[Nscalar];
  vecType  vecInput[Nvec], vecOutp[Nvec], vecZeroes= vecType( 0UL ) ;  
  static Timer<nanoseconds> timer;

  scalarLong* resVector= new scalarLong[nRepetitions];
  scalarLong* resScalar= new scalarLong[nRepetitions];
  
  for( int i= 0; i< nRepetitions; i++) {
    resVector[i]= 0;
    resScalar[i]= 0;
  }

  // Prepare scalar array
  for( int i= 0; i< Nscalar; i++){
    arrA[i]= i;
    arrOutVector[i]= 0UL;
    arrOutScalar[i]= 0UL;
  }

  for( int ind= 0; ind< Nvec; ind++){
    vecOutp[ind] = vecZeroes;
  }
  
  vecType  vecState[Nvec];
  PrepareInputVecArray( vecInput, arrA ); // 3rd=verbose ( def = false )
  //******************
  
  constexpr int lastInd= (Nscalar-1) % VecWidth;
  // int   slot = Nscalar; // For use in changing the values ...

  bool varyInput = true;

  // --- Scalar  Version ---  Start -------------------  
  timer.Start();
  for( int j= 0; j < nRepetitions; j++)  
  {
     // for( int i= 0; i< Nscalar; i++)  arrA[i]= i;
     int indChanged = 0;
     if( varyInput ) {
        indChanged = j % Nscalar;
        arrA[ indChanged ] = 0;
     }

     if( enableModulo ) {
        SimplePartialSumWithModulo( arrA, Nscalar, arrOutScalar );
     } else {
        SimplePartialSum( arrA, Nscalar, arrOutScalar );
     }
     
     resScalar[j]= arrOutScalar[Nscalar-1];
     if( varyInput ) {     
        arrA[ indChanged ] = indChanged;
     }
  }
  double scalarTime = timer.Elapsed();  
  // --- Scalar  Version ---  Start -------------------  
  
  // --- Vector  Version ---  Start -------------------  
  timer.Start();
  for( int j= 0; j< nRepetitions; j++)  
  {
     // slot--; if( slot == 0 ) { slot = Nscalar; }
     int iChanged = j % Nscalar;

     int indV = iChanged >> 2;
     int slot = iChanged & 3;
     vecCore::Set( vecInput[indV], slot, 0UL );
        
     PartialSumVec<vecType, Nvec, enableModulo>( vecInput, vecOutp );
     //*********** ---------------------------
     resVector[j]= vecCore::Get( vecOutp[Nvec-1], lastInd );

     if( varyInput ) {
        vecCore::Set( vecInput[indV], slot, iChanged );
     }
  }
  double vecTime = timer.Elapsed();

  ConvertVecToArray(  vecOutp, arrOutVector, Nscalar, true );  // last is 'verbose'
  // --- Vector  Version ---   End  -------------------
  
  // std::cout << " Results for N (scalar) = " << Nscalar
  //          << " N-vector = " << Nvec << std::endl;
  printf( " Results for N= %4d ( vector len= %3d ) - number of repetitions=  ",
          Nscalar, Nvec);
  if( nRepetitions > 100000 ||
     nRepetitions == ( nRepetitions / 1000) * 1000 ) {
     printf( "%8d K \n", nRepetitions / 1000);
  } else {
     printf( "%6d \n", nRepetitions );
  }     

  bool verbose= false;
  if( verbose ) { 
     // Print last result .. 
     ReportSumResults( arrA, arrOutScalar, vecOutp,  // arrOutVector
                       Nscalar );
  }
  
  // Check last result
  int nBad= 0;
  
  for( int i= 0; i< Nscalar; i++)
  {
    if( arrOutVector[i] != arrOutScalar[i] ) { 
      nBad ++; 
      std::cerr << " Last result> Error in location " << i 
                << "  value = " << arrOutVector[i]
                << "  expected = " << arrOutScalar[i] << std::endl;
      // fprintf( stderr, " Error in location %d : value= %ld  - expected %d \n",
      //                 i, arrA[i], i+1 );
    }
  }

  if( nBad > 0 ) {
     std::cerr << " Problem result found in comparison with simple sequential method. " << std::endl;
     std::cerr << " -- Summary of comparison: Good = " << Nscalar - nBad
               << " Bad = " << nBad << " out of " << Nscalar << " . " << std::endl;
     std::cerr << "    Percentage bad = " << 100.0 * nBad / Nscalar << " % "
               << std::endl;
  }

  // Check all sums
  int nBad2= 0;
  for( int i= 0; i< nRepetitions; i++)
  {
    if( resVector[i] != resScalar[i] ) { 
      nBad2 ++; 
      // std::cerr << " Error in location " << i 
      //       << "  value = " << arrA[i] 
      //       << "  expected = " << i+1 << std::endl;
      fprintf( stderr, " Error in sum of repetition  %d : value= %llu  - expected %llu \n",
                       i, resVector[i], resScalar[i] );
    }
  }

//  std::cout << "Time of templated method= " << std::setprecision(3) << vecTime / nRepetitions << " ns "
//            << " (per trial ) " << std::endl
//            << "     of simple    method= " << scalarTime / nRepetitions << " ns " << std::endl;

  fprintf( stdout, "Time of  vector  method= %5.1f ns (per trial)\n", vecTime / nRepetitions);
  fprintf( stdout, "     of  scalar  method= %5.1f ns (per trial)\n", scalarTime / nRepetitions);

  delete[] resVector;
  delete[] resScalar; 
  
  return nBad + nBad2;
}
#endif

void ReportSumResults( const scalarLong arrInput[],
                       const scalarLong outScalar[],  // Reference result
                       const vecType    vecSum[Nvec], // New (vec) result
                       int   lenArr )
{
  vecType tokenVT;
  constexpr int vecWidth= vecCore::VectorSize<decltype(tokenVT)>();
  
  printf( "   %3s    %10s   %10s  %10s \n", "I", "Input", "Scalar-psum", "Vec-par/sum");
  printf( "========================================================================\n");
  for( int i= 0; i< lenArr; i++)
  {
    scalarLong vecRes= vecCore::Get( vecSum[i/vecWidth], i % vecWidth );
    // std::cout << " [ " << i << " scalar value  = " << arrA[i];
    printf( " [ %3d ]  %10llu   %10llu  %10llu", i, arrInput[i], outScalar[i], vecRes);
    if( outScalar[i] != vecRes ) {
       printf ( " Difference = %lld ", (long long) outScalar[i] - vecRes);
    }
    printf( "\n" );
  }
  printf( "========================================================================\n"); 
}

int main( int argc, char **argv)
{
  // if( argc == 1) N = atoi( argv[0] );
  TestPartialSum();

  TestPartialSum( true ); // i.e. enable Modulo

#ifdef  ENABLE_BENCHMARK
  int nRepetitions= 100000;

  std::cout << " Argc = " << argc;
  if( argc > 1 ) {
     int nrepsIn = atoi( argv[1] );

     std::cout << " numRepetitions = Argv[1] = " << nrepsIn ;
     if( nrepsIn > 0 )
        nRepetitions = nrepsIn;
  }

  BenchmarkPartialSums(nRepetitions);
#endif
}
