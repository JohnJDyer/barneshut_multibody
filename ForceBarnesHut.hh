#ifndef _FORCEBARNESHUT_HH_
#define _FORCEBARNESHUT_HH_

#include "ForceCalculator.hh"
#include "MortonKeyCalculator.hh"
#include "Body.hh"
#include <vector>
#include <memory>
#include <tuple>
#include <iostream>
#include <array>
#include <utility>

#define BISECTION   0
#define LEFT_IS_LEAF  1
#define RIGHT_IS_LEAF 2

struct I_Body_2 : public Body {
  //using Body::Body;
  //I_Body_2(double x, double y, double vx, double vy, double m):Body(x,y,vx,vy,m){}

  I_Body_2(
    std::array<I_Body_2*, 4>&,
    std::array<Body*,   4>&,
    double,
    double,
    double,
    double);

  ~I_Body_2(){
    for(auto i = 0; i < 4; i ++ ){
    delete _i_subbodies[i];
    }
  };

  std::array<I_Body_2*, 4> _i_subbodies;
  std::array<Body*,   4>   _subbodies;
  double diameter;

  void print_tree(int n);
  void apply_force(Body*, double);
};


struct Compare {
  Compare(int level, MortonKeyCalculator *mkcalc): level(level), mkcalc(mkcalc) {}

  bool operator()(Body &a, int i ){
    int bits = ((mkcalc->y(a.y()) & level) >0) * 2 + ((mkcalc->x(a.x()) & level) >0) * 1;
    return bits < i;
  }
  int level;
  MortonKeyCalculator* mkcalc;

};

class ForceBarnesHut : public ForceCalculator {

 public:
  ForceBarnesHut(Body *body, int N, double theta);
  ForceBarnesHut() = delete;

  virtual void operator() (Body *pulled);
  virtual ~ForceBarnesHut(){delete tree;};

  I_Body_2* construct_tree(Body*, Body*, uint32_t, double);
  void print_tree();

 private:
  Body *bodies_;
  const int       N_;
  double        theta_;

  MortonKeyCalculator mkcalc;

  I_Body_2* tree;
};

#endif
