
//Scalar Ln(n!) function for Discrete Poisson PDF
template <typename DerivedT>
template <typename BackendT>
VECCORE_ATT_HOST_DEVICE typename BackendT::Double_v
VecRNG<DerivedT>::LogFactorialScalar(typename BackendT::Double_v LG,
                         typename BackendT::Double_v N,
                          typename BackendT::Double_v z){

   using Double_v = typename BackendT::Double_v;

   Double_v lct[9] = {
    //  1.000000000190015,
    // 76.18009172947146,
  //  -86.50532032941677,
  //   24.01409824083091,
  //   -1.231739572450155,
  //    0.1208650973866179e-2,
      	0.99999999999980993227684700473478,
       	676.520368121885098567009190444019,
       	-1259.13921672240287047156078755283,
       	771.3234287776530788486528258894,
       	-176.61502916214059906584551354,
       	12.507343278686904814458936853,
       	-0.13857109526572011689554707,
       	9.984369578019570859563e-6,
       	1.50563273514931155834e-7};

 Double_v ln_sqrt_2_pi = math::Log(math::Sqrt(2*M_PI));
if (z < 0.5) {

  Double_v p=(M_PI/math::Sin(M_PI*z));
  Double_v w=LogFactorialScalar<BackendT>(LG,N,1.0-z);
return p-w;
}else{
z = z-1.0;
 // Base of the Lanczos exponential
Double_v base =z+LG+0.5;
Double_v sum = 0;
// We start with the terms that have the smallest coefficients and largest
// denominator.
for(int i=6; i>=1; i--) {
sum += lct[i] / (z + ( (Double_v) i) );
}
sum += lct[0];
Double_v res=(ln_sqrt_2_pi + math::Log(sum))-base;
Double_v t=math::Log(base)*(z+0.5);
return res+t;
}
}

//Vectorized Ln(n!) function for Discrete Poisson PDF
template <typename DerivedT>
template <typename BackendT>
VECCORE_ATT_HOST_DEVICE typename BackendT::Double_v
VecRNG<DerivedT>::LogFactorial(typename BackendT::Double_v LG,
                        typename BackendT::Double_v N,
                          typename BackendT::Double_v z){
  using Double_v = typename BackendT::Double_v;
 int vsize = VectorSize<Double_v>();

    Double_v lct[9] = {
     //  1.000000000190015,
     // 76.18009172947146,
   //  -86.50532032941677,
   //   24.01409824083091,
   //   -1.231739572450155,
   //    0.1208650973866179e-2,
       	0.99999999999980993227684700473478,
        	676.520368121885098567009190444019,
        	-1259.13921672240287047156078755283,
        	771.3234287776530788486528258894,
        	-176.61502916214059906584551354,
        	12.507343278686904814458936853,
        	-0.13857109526572011689554707,
        	9.984369578019570859563e-6,
        	1.50563273514931155834e-7 };


 Double_v ln_sqrt_2_pi = math::Log(math::Sqrt(2*M_PI));
 Double_v p;
 Double_v w;

 Mask<Double_v> D (z < Double_v(0.5));
 if(MaskFull(D)) {
    p=(M_PI/math::Sin(M_PI*z));
    w=LogFactorial<BackendT>(LG,N,1.0-z);
   return p-w;
 }else{
 if(MaskEmpty(D)) {
 z = z-1.0;
 Double_v base =z+LG+0.5;
 Double_v sum = 0;

 for(int i=6; i>=1; i--) {
 sum += lct[i] / (z + ( (Double_v) i) );
 }
 sum += lct[0];
 Double_v res=(ln_sqrt_2_pi + math::Log(sum))-base;
 Double_v t=math::Log(base)*(z+0.5);
 return res+t;
}else{
  Double_v r;
for(int i = 0 ; i < vsize ; ++i) {
  r[i]=LogFactorialScalar<ScalarBackend>(LG[i],N[i],z[i]);
}
return r;
}
}
}
