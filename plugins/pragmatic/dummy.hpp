#ifndef DUMMY_HPP
#define DUMMY_HPP

#include <chrono>
#include <thread>

double dummy_call (int x, double y)
{
  double volatile dummy;

  for (size_t i = 0; i < 50; ++i)
  {
    for (size_t j = 0; j < 50; ++j)
    {
      dummy = atan(pow(x, i)) + sin(exp(y)) / (j*sqrt(y));
    }  
  }

  return dummy;
}
#endif
