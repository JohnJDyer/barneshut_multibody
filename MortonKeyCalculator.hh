#ifndef __MORTONKEYCALCULATOR_H__
#define __MORTONKEYCALCULATOR_H__

#include <ostream>
#include "Body.hh"
#include <tuple>

class MortonKeyCalculator{
  public:

    MortonKeyCalculator(Body *bodies, int N);

    uint32_t x(double xval) const{
      return (uint32_t)(0xffffffff * (xval - xmin) / max_range);
    };

    uint32_t y(double yval) const{
      return (uint32_t)(0xffffffff * (yval - ymin) / max_range);
    }

    double cellWidth(int level) const{
      return max_range / (1 << (level -1));
    }

    std::ostream& printKey(std::ostream &os, const Body &b) const;

    bool operator() (const Body &a, const Body &b) const;

//    int radix_dist(int, int);
//    std::tuple<uint32_t, uint32_t> radix_prefix(int a_idx, uint32_t length);


private:
    double xmin;
    double ymin;
    double max_range;
    Body  *bodies;
    int    N_bodies;
};

#endif
