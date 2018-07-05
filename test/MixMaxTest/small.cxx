//  Versions for use with compiler explorer
//   e.g. at godbolt.org



#include <cstdint>     // For int64_t and similar 
#include <algorithm>   // For std::min

void SimplePartialSum(const uInt64_t inArr[], int num, uInt64_t outArr[], uInt64_t overflowArr[] = 0 )
{
  outArr[0] = inArr[0];
  uInt64_t  overflow= 0;
  for( int i= 1; i< num; i++){
    outArr[i] = outArr[i-1] + inArr[i];
    if( outArr[i] < outArr[i-1] )
       overflow++;
    if( overflowArr )
       overflowArr[i] = overflow;
  }
}


void
PartialSumV2( const uint64_t inpArr[], uint64_t outArr[], unsigned int N )
{
  for( unsigned int i=1; i< N; i += i )
  {
     outArr[i]= inpArr[i];
  }
  for( unsigned int i=1; i< N; i += i )
  {
    unsigned int len= i; 
    unsigned int skip= i+i;
    for( unsigned int j=i; j< N; j += skip )
    {      
       uint64_t prev= outArr[j-1];
       unsigned int kTop = std::min( j + len, N );
       for( unsigned int k=j; k< kTop; k++ )
       {
          outArr[k] += prev;
       }
    }
  }
}

// ==============================

template< typename dataType, unsigned int N >
void
PartialSumV2( const dataType inpArr[], dataType outArr[] )
{
  for( unsigned int i=1; i< N; i += i )
  {
     outArr[i]= inpArr[i];
  }
  for( unsigned int i=1; i< N; i += i )
  {
    unsigned int len= i; 
    unsigned int skip= i+i;
    for( unsigned int j=i; j< N; j += skip )
    {      
       dataType prev= outArr[j-1];
       unsigned int kTop = Min( j + len, N );
       for( unsigned int k=j; k< kTop; k++ )
       {
          outArr[k] += prev;
       }
    }
  }
}
