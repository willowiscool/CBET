#ifndef BEAM_H
#define BEAM_H

#include <vector>
#include <tuple>
#include <cstddef>
#include <omp.h>

using namespace std;

struct Crossing {
	double x;
	double z;
	size_t boxesx;
	size_t boxesz;
	double area_ratio;
	double dkx;
	double dkz;
	double dkmag;
	double i_b;
	double w_mult;
};

struct Ray {
	vector<Crossing> crossings;
	double x0;
	double z0;
	double cx0;
	double cz0;
	double kx0;
	double kz0;
};

struct Beam {
	vector<Ray> rays;
	vector<tuple<omp_lock_t, vector<tuple<size_t, size_t>>>> marked;
	vector<tuple<bool, tuple<size_t, size_t>>> raystore;
};

Beam beam1();
Beam beam2();

#endif
