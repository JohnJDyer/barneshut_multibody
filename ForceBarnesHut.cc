#include "ForceBarnesHut.hh"
#include "Timer.hh"

#ifdef _OPENMP
#include <parallel/algorithm>
#include <parallel/numeric>
#else
#include <algorithm>
#endif

#include <iostream>
#include <iterator>     // std::ostream_iterator
#include <vector>       // std::vector
#include <algorithm>    // std::copy
#include <memory>
#include <tuple>
#include <utility>
#include <array>

//#include "MortonKeyCalculator.hh"
#include "MortonKeyCalculator.cc"

using namespace std;


I_Body_2::I_Body_2(array<I_Body_2*, 4>& i_subbodies, array<Body*, 4>& subbodies, double diameter, double x, double y, double m):
	Body(x,y,0.0,0.0,m), _i_subbodies(i_subbodies), _subbodies(subbodies), diameter(diameter)
{

};




I_Body_2* ForceBarnesHut::construct_tree(Body *left, Body *right, uint32_t level, double diameter) {
	//cout << level << " " << right -left << endl <<  flush;

	array<Body*, 5> bounds 					= {} ;
	array<Body*, 4> subbodies 			= {};
	array<I_Body_2*, 4> i_subbodies = {};

	bounds[0] = left;
	bounds[4] = right;

//	auto compare = [&](Body &a, int i) -> bool {
//			int bits = ((mkcalc.y(a.y()) & level) >0) * 2 +
//							   ((mkcalc.x(a.x()) & level) >0) * 1;
//			return bits < i;
//	};

	auto compare = Compare(level, &mkcalc);

	# pragma omp parallel for
  for(int i = 1; i < 4; i++){
		bounds[i] = lower_bound(left, right, i, compare); //push this down.

		//cout << "\t" << x - left << endl;
	}

	I_Body_2 *t_0 = NULL, *t_1 = NULL, *t_2 = NULL, *t_3 = NULL;

	auto sub_level    = level >> 1;
	auto sub_diameter = diameter / 2.0;

	int i = 0;
	const int make_child = 1;

	//const int parallel_thresh = 10; //"time_setup": 0.986188
	//const int parallel_thresh = 15; //"time_setup": 0.847496,
	//const int parallel_thresh = 20; //"time_setup": 0.791603
	const int parallel_thresh = 25; //"time_setup": 0.592091,
	//const int parallel_thresh = 30; //"time_setup": 0.703603,
	//const int parallel_thresh = 40; //"time_setup": 1.10037



	if(bounds[i+1] - bounds[i] > make_child){
		#pragma omp task shared(t_0) if(bounds[i+1] - bounds[i] > parallel_thresh)
		t_0 = construct_tree(bounds[i], bounds[i+1], sub_level, sub_diameter);
	}

	i = 1;
	if(bounds[i+1] - bounds[i] > make_child){
		#pragma omp task shared(t_1) if(bounds[i+1] - bounds[i] > parallel_thresh)
		t_1 = construct_tree(bounds[i], bounds[i+1], sub_level, sub_diameter);
	}

	i = 2;
	if(bounds[i+1] - bounds[i] > make_child){
		#pragma omp task shared(t_2) if(bounds[i+1] - bounds[i] > parallel_thresh)
		t_2 = construct_tree(bounds[i], bounds[i+1], sub_level, sub_diameter);
	}

	i = 3;
	if(bounds[i+1] - bounds[i] > make_child){
		#pragma omp task shared(t_3) if(bounds[i+1] - bounds[i] > parallel_thresh)
		t_3 = construct_tree(bounds[i], bounds[i+1], sub_level, sub_diameter);
	}

	for(int i = 0; i < 4; i++){
		if(bounds[i+1] - bounds[i] == 1){ //if fewer actually
			subbodies[i] = bounds[i];
		} else {
			subbodies[i] = NULL;
		}
	}

	#pragma omp taskwait

	i_subbodies[0] = t_0;
	i_subbodies[1] = t_1;
	i_subbodies[2] = t_2;
	i_subbodies[3] = t_3;


	double t_m=0.0, t_x=0.0, t_y=0.0;

	for(int i=0; i < 4; i ++){
		if (i_subbodies[i] != NULL){
			auto & sub = i_subbodies[i];
			double l_m = sub-> m();
			t_x += l_m * sub -> x();
			t_y += l_m * sub -> y();
			t_m += l_m;
		} else if (subbodies[i] != NULL) {
			auto & sub = subbodies[i];
			double l_m = sub-> m();
			t_x += l_m * sub -> x();
			t_y += l_m * sub -> y();
			t_m += l_m;
		}
	}

	t_x /= t_m;
	t_y /= t_m;

	return new I_Body_2(i_subbodies, subbodies, diameter, t_x, t_y, t_m); // this can be moved up

}


void I_Body_2::apply_force(Body *pulled, double theta) {

	auto dx = pulled->x() - this->x();
	auto dy = pulled->y() - this->y();

	auto dist_squared     = dx*dx + dy*dy;
	auto diam_squared     = this->diameter * this->diameter;

	if(diam_squared / dist_squared  < theta * theta){
		pulled->accGravityFrom(*this);
	} else {
		for(int i = 0; i < 4; i++){
			if(_i_subbodies[i] != NULL){
				_i_subbodies[i] -> apply_force(pulled, theta);
			} else if(_subbodies[i] != NULL) {
				pulled->accGravityFrom(*_subbodies[i]);
			}
		}
	}

};

void I_Body_2::print_tree(int n) {

	Body* test = (Body*) this;

	for(int i = 0; i < n; i++){
		cout << "\t";
	}

	cout << "i:x,y,m:" << test ->x() << "\t"<< test ->y () << "\t" << test ->m() << "\n";

	for(int i = 0; i < 4; i++){
		if(_i_subbodies[i] != NULL){
			_i_subbodies[i] -> print_tree(n+1);
		} else if(_subbodies[i] != NULL) {
			for(int a = 0; a < n+1; a++){
				cout << "\t";
			}

			cout << "l:x,y,m:" << _subbodies[i]->x() << "\t"<< _subbodies[i]->y() << "\t" << _subbodies[i]->m() << "\n";

		}
	}

};



//Build the quad tree and do any other needed initialization here.
//Quad tree construction should be parallelized
ForceBarnesHut::ForceBarnesHut(Body *bodies, int N, double theta):
	bodies_(bodies), N_(N), theta_(theta), mkcalc(bodies_, N)
{

	__gnu_parallel::sort(bodies, bodies + N, mkcalc);

	#pragma omp parallel
	#pragma omp single
	{
		//cout << "start diameter: " <<  mkcalc.cellWidth(1) <<  "\t" << (1u << 31u) << endl;
		tree =construct_tree(bodies, bodies+N, 1u << 31u, mkcalc.cellWidth(1));
	};

	//tree ->print_tree(0);

};


void ForceBarnesHut::operator()(Body *pulled){

	tree -> apply_force(pulled, theta_);

	//Replace the naive code below with an efficient implementation
	//of Barnes-Hut.
	//	for(int i = 0; i < N_; i++){
	//		pulled->accGravityFrom(bodies_[i]);
	//	}

};


//	for(auto & b = left; left != right; left++){
//		mkcalc.printKey(cout, *b);
//		cout << " ";
//		cout << compare(*b,0);
//		cout << compare(*b,1);
//		cout << compare(*b,2);
//		cout << compare(*b,3);
//		cout << endl;
//		//cout << compare(*b,i);
//	}
