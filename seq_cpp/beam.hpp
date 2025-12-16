#ifndef BEAM_H
#define BEAM_H

#include <vector>
#include <tuple>
#include <cstddef>

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
	std::vector<Crossing> crossings;
	double x0;
	double z0;
	double cx0;
	double cz0;
	double kx0;
	double kz0;
};

struct Beam {
	std::vector<Ray> rays;
	std::vector<std::vector<std::tuple<size_t, size_t>>> marked;
	std::vector<std::tuple<bool, std::tuple<size_t, size_t>>> raystore;
};

Beam beam1();
Beam beam2();

#endif
