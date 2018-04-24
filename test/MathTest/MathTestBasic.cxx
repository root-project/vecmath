#include <VecCoreLib/Math.h>

#include <cmath>
#include <cstdlib>
#include <iostream>

int main() {
  std::cout << "Exp:" << std::endl;
  std::cout << vecMath::FastExp(10.0) << std::endl;
  std::cout << vecMath::FastExp((vecCore::backend::VcVector::Double_v)10.0)
            << std::endl;

  std::cout << "Log:" << std::endl;
  std::cout << vecMath::FastLog(10.0) << std::endl;
  std::cout << vecMath::FastLog((vecCore::backend::VcVector::Double_v)10.0)
            << std::endl;

  std::cout << "Pow:" << std::endl;
  std::cout << vecMath::FastPow(12.0, 0.3) << std::endl;
  std::cout << vecMath::FastPow((vecCore::backend::VcVector::Double_v)12.0,
                                (vecCore::backend::VcVector::Double_v)0.3)
            << std::endl;
  std::cout << "std::pow " << std::pow(12.0, 0.3) << std::endl;

  srand(time(NULL));
  double x = rand();
  std::cout << "IntPow:" << std::endl;
  std::cout << vecMath::IntPow(x, 3) << std::endl;
  std::cout << vecMath::IntPow((vecCore::backend::VcVector::Double_v)x, 3)
            << std::endl;
  std::cout << "std::pow " << std::pow(x, 3) << std::endl;
  return 0;
}
