#include <iostream>
#include "ProgPow.h"

// g++ test.cpp ProgPow.cpp; ./a.out
int main(void) {
   std::cout << ProgPow::getKern(123, ProgPow::KERNEL_CUDA);
}
