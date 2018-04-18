#include <VecCoreLib/Math.h>

#include <iostream>

int main(){
  std::cout << "Exp:" << std::endl;
  std::cout << vecMath::FastExp(10.0) << std::endl;
  std::cout << vecMath::FastExp((vecCore::backend::VcVector::Double_v) 10.0) << std::endl;

  std::cout << "Log:" << std::endl;
  std::cout << vecMath::FastLog(10.0) << std::endl;
  std::cout << vecMath::FastLog((vecCore::backend::VcVector::Double_v) 10.0) << std::endl;
  return 0;
}
