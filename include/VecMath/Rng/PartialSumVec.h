#ifndef  PartialSumVec_h
#define  PartialSumVec_h

// #define  __AVX2__   1  // Instruct UME::SIMD to use AVX2 ... !?

#include "VecCore/VecCore"

// #include "VecCore/Backend/VcSimdArray.h"
// #include "VecCore/VecCore"
// #include "VecCore/Backend/UMESimdArray.h"

//  Constants - until moved to common header
//
constexpr uint64_t gMers61 = 0x1fffffffffffffff;

constexpr int VecWidth=4;  // For AVX ... 64 x 4 = 256 bits
// static constexpr int vecSize= VectorSize<decltype(dbVec)>();

#if 0
template< typename dataType>
std::ostream& operator<< ( std::ostream& os,
                           const dataType & vecVal )
// Streaming operator
{
  constexpr int vecSize= vecCore::VectorSize<decltype(vecVal)>();
  for( int i=0 ; i< vecSize ; i++ )
     os // << " "
        << const_cast<typename vecCore::TypeTraits<dataType>::UInt64_t>(vecCore::Get( vecVal, i ))
        << " "
        ;

  return os;
}
#else
#include <iomanip>

std::ostream& operator<< ( std::ostream& os,
                           const vecCore::backend::UMESimdArray<4>::UInt64_v & vecVal )
// Streaming operator
{
  constexpr int vecSize= vecCore::VectorSize<decltype(vecVal)>();
  for( int i=0 ; i< vecSize ; i++ )
     os << " " << std::setw(3) << vecCore::Get( vecVal, i );

  return os;
}
#endif

// ---------------------------------------------------------------
// If two numbers are both 0 <= a, b < P
//   the following method can be used for their sum Mod p
//
template< typename dataType>
dataType 
SumModuloP( const dataType a, const dataType b, const uint64_t P )
{
  vecCore::Mask<dataType> overP= a > (dataType) P;
  MaskedAssign( a, overP, a-P );
}

// ---------------------------------------------------------------
//
//

template< typename dataType, bool largeVal = false >
// inline
VECCORE_FORCE_INLINE
dataType 
FastModuloMers61( const dataType a )
{
  auto mask61 = dataType( gMers61 );
  
  auto low = a & mask61; // dataType( gMers61 );  //  Max = gMers61 = 0x0fff ffff ffff ffff
  auto hi  = a >> 61;      //  Max = 7 
  auto s   = low + hi;     //  Max = 0x1000 0000 0000 0006

  auto low2 = s & mask61; // dataType( gMers61 );
  auto hi2  = s >> 61;     //  Max = 1
  return low2 + hi2;       
}

template< typename dataType>
void PrintArrLine( const char* label, const dataType inpArr[], int N );

using std::cout;
using std::endl;

// #define VERBOSE  1

// ---------------------------------------------------------------
//
//
template< typename dataType, unsigned int N, bool enableMod >
void
PartialSumVec( const dataType inpArr[N], dataType outArr[N] )
{
   using vecCore::Set;
   using vecCore::Get;
   // using ScalarType = uint64_t; // or Scalar<dataType> ...
   using ScalarType = typename vecCore::TypeTraits<dataType>::ScalarType;   

   using Mask_v = typename vecCore::TypeTraits<dataType>::MaskType;

   const Mask_v MaskStep2 // = { false, false, true, true };
                          ( false, false, true, true );
   constexpr int vecSize= vecCore::VectorSize<dataType>(); // decltype(dummy)>();

#ifdef VERBOSE
   dataType sum1keep[N];
   constexpr bool verbose= true;
   dataType shift1keep[N], shift2keep[N];
#else   
   constexpr bool verbose= false;
#endif

   // using Mask_v = Vc::SimdMaskArray<ScalarType, VecWidth>;  // Works for Vc - not general
   // using Mask_v = vecCore::SimdMaskArray<ScalarType, VecWidth>;  // Failes for all
   
   if( verbose ) {
     cout << " vecSize = " << vecSize << endl;
     PrintArrLine( "Array In", inpArr, N );
   }

   for( int i= 0; i<N; ++i) {
      // Shift one - from even to odd only      
      dataType  shift1( 0UL ), sum1( 0UL );
      dataType tshift2( 0UL );
      
      Set( shift1, 1, Get( inpArr[i], 0 ));
      Set( shift1, 3, Get( inpArr[i], 2 ));
      sum1 = inpArr[i] + shift1;

      ScalarType  val1= Get( sum1, 1 );
      tshift2= vecCore::Blend( MaskStep2, dataType(val1), dataType( 0L ) );

      outArr[i] = sum1 + tshift2;
      // auto res1 = sum1 + tshift2;    // Idea 2.5
      // Sum of up to four values

      // outArr[i] = sum1[i] + tshift2;    // Gives interesting compiler errors
      
#ifdef VERBOSE
      shift1keep[i] = shift1;
      sum1keep[i] = sum1; 
      shift2keep[i] = tshift2;
      // cout << tshift2;
#endif

// #define SEPARATE_LOOPS 1
      
#ifdef SEPARATE_LOOPS
   }
   
#ifdef VERBOSE
   if( verbose ) {
      PrintArrLine( "Shift1   ", shift1keep, N ); /*<dataType>*/
      PrintArrLine( "Sum1     ", sum1keep,   N );
      PrintArrLine( "Shift2   ", shift2keep, N );
      PrintArrLine( "Sum2     ", outArr, N );
   }
#endif   
   // Expected result (til now): 
   // Global: [ 0   1    2     3   ]    [  4   5    6     7   ]
   //  'i'     ------- 0 ----------    ---------- 1 ----------
   // Local:  [ 0   1    2     3   ]    [  0   1    2     3   ]  
   // Sums:   [ 0  0+1   0-2  0-3  ]    [  4  4+5  4-6   4-7  ]

   // Global: [ 8   9   10    11   ]    [ 12   13    14    15   ]
   //  'i'     ------- 2 ----------    ---------- 3 ----------
   // Local:  [ 0   1    2     3   ]    [  0    1     2     3   ]  
   // Sums:   [ 8  8+9  8-10  8-11 ]    [ 12  12-13 12-14 12-15  ]

   // "Version 2" - does all the remaining work, but requires
   //    results of previous iteration in next one.
   for( int i= 1; i<N; i++)
   {
      // If loops are 'separated' the following works:      
      ScalarType  sumPrevious= Get( outArr[i-1], 3 );      
#else
      ScalarType  sumPrevious= i > 0 ? Get( outArr[i-1], 3 )  : 0UL ;
#endif
      auto res = outArr[i] + dataType( sumPrevious );
      if( enableMod ){
         outArr[i] = FastModuloMers61( res );
      } else {
         outArr[i] = res;
      }
   }
}

// ---------------------------------------------------------------
template< typename dataType>
void PrintArrLine( const char* label, const dataType arrayVec[], int N )
{
   using std::cout;
   // cout << endl;
   cout << std::setw(9) << label << " = " ;
   for( int i= 0; i<N; ++i)
      cout << arrayVec[i]; // << " - ";
   cout << std::endl;
}

/***
template< typename dataType, unsigned int N >
void
PartialSumV3( dataType a[] ) // , unsigned int N )  // [N] // InOut
{
  int intermed= 4;
  for( unsigned int i=1; i< intermed; i += i )
  {
    unsigned int len= i; 
    unsigned int skip= i+i;
    for( unsigned int j=i; j< N; j += skip )
    {      
       dataType prev= a[j-1];
       unsigned int kTop = Min( j + len, N );
       for( unsigned int k=j; k< kTop; k++ )
       {
          a[k] += prev;
       }
    }
  }
  for( unsigned int i=intermed; i< N; i += i )
  {
    unsigned int len= i; 
    unsigned int skip= i+i;
    for( unsigned int j=i; j< N; j += skip )
    {      
       dataType prev= a[j-1];
       unsigned int kTop = Min( j + len, N );
       for( unsigned int k=j; k< kTop; k++ )
       // Must be fully vectorised (modulo tail ? )
       {
          a[k] += prev;
       }
    }
  }  
}
***/

template< typename dataType, unsigned int N>
void
AltPartialSumModuloP( dataType a[], const int P ) // , unsigned int N )  // [N] // InOut
{
  using Mask_v = typename vecCore::TypeTraits<dataType>::MaskType;
  for( unsigned int i=1; i< N; i += i )
  {
    unsigned int len= i; 
    unsigned int skip= i+i;
    for( unsigned int j=i; j< N; j += skip )
    {      
       dataType prev= a[j-1];
       unsigned int kTop = vecCore::math::Min( j + len, N );
       for( unsigned int k=j; k< kTop; k++ )
       {
          a[k] += prev;
          Mask_v overP= a[k] > P;
          MaskedAssign( a[k], overP, a[k]-P );
          //  or use a FastModulo(x,P61) ??
       }
    }
  }
}

#endif
