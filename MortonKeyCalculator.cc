#ifdef _OPENMP
#include <parallel/algorithm>
#include <parallel/numeric>
#else
#include <algorithm>
#endif

#include "MortonKeyCalculator.hh"
#include<iostream>
#include<climits>
#include<tuple>
#include<utility>
#include <cmath>
#include <limits>
#include <omp.h>

bool x_comp(Body &a, Body &b){
  return a.x() < b.x();
}

bool y_comp(Body &a, Body &b){
  return a.y() < b.y();
}

MortonKeyCalculator::MortonKeyCalculator(Body *bodies, int N){
  this -> bodies   = bodies;
  this -> N_bodies = N;

  auto max_x = __gnu_parallel::max_element(bodies, bodies +N , x_comp);
  auto max_y = __gnu_parallel::max_element(bodies, bodies +N , y_comp);
  auto min_x = __gnu_parallel::min_element(bodies, bodies +N , x_comp);
  auto min_y = __gnu_parallel::min_element(bodies, bodies +N , y_comp);

  this -> xmin = min_x -> x();
  this -> ymin = min_y -> y();

  this -> max_range = std::max(
    max_y -> y() -  min_y -> y(),
    max_x -> x() -  min_x -> x());
};

std::ostream& MortonKeyCalculator::printKey(std::ostream &os, const Body &b) const{
  uint32_t bx = x(b.x());
  uint32_t by = y(b.y());

  for(int i = 31; i >= 0; i--){
    os<<((by >> i) % 2);
    os<<((bx >> i) % 2);
  }

  os<<std::endl;
  return os;
}

inline bool less_msb(uint32_t x, uint32_t y) {
  return (x < y) and (x < (x ^ y));
}

bool MortonKeyCalculator::operator() (const Body &a, const Body &b) const{
  uint32_t _a[2] = {y(a.y()), x(a.x())};
  uint32_t _b[2] = {y(b.y()), x(b.x())};

  int      j = 0;
  uint32_t x = 0;

  for (int k = 0; k < 2; k++) {
    uint32_t y = _a[k] ^ _b[k];
    if (less_msb(x, y)) {
      j = k;
      x = y;
    }
  }

  return _a[j] < _b[j];

}

//int MortonKeyCalculator::radix_dist(int a_idx, int b_idx){
//
//  if( a_idx < 0 or a_idx > N_bodies-1 or b_idx < 0 or b_idx > N_bodies-1)
//    return -1;
//
//  auto a = bodies[a_idx];
//  auto b = bodies[b_idx];
//
//  uint32_t _a[2] = {y(a.y()), x(a.x())};
//  uint32_t _b[2] = {y(b.y()), x(b.x())};
//
//  auto x_count = __builtin_clz(_a[1] ^ _b[1]);
//  auto y_count = __builtin_clz(_a[0] ^ _b[0]);
//
//  return std::min(x_count, y_count) * 2 + (y_count > x_count);
//
//}
//
//std::tuple<uint32_t, uint32_t>  MortonKeyCalculator::radix_prefix(int idx, uint32_t length){
//
//  auto a = bodies[idx];
//
//  uint64_t __y = y(a.y());
//  uint64_t __x = x(a.x());
//
//  return std::make_tuple(
//    __y >> (32 -(length >> 1) - (length % 2 == 1)),
//    __x >> (32 -(length >> 1)                    )
//  );
//
//}
