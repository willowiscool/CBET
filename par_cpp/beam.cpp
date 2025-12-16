#include <cstddef>
#include <omp.h>
#include "beam.hpp"
#include "consts.hpp"

using namespace std;

Beam beam1() {
	double dz = (consts::BEAM_MAX_Z-consts::BEAM_MIN_Z)/((double)consts::NRAYS-1.0);
	vector<Ray> rays;
	rays.reserve(consts::NRAYS);
	for (size_t i = 0; i < consts::NRAYS; i++) {
		double z0 = (double)i * dz + consts::BEAM_MIN_Z + consts::OFFSET1;
		rays.push_back(Ray {
			vector<Crossing>(), // crossings
			consts::XMIN, // x0
			z0,
			consts::XMIN, // cx0
			z0 + consts::CHILD_OFFSET, // cz0
			consts::DIR1[0], // kx0
			consts::DIR1[1], // kz0
		});
	}
	vector<tuple<omp_lock_t, vector<tuple<size_t, size_t>>>> marked
		(consts::NX*consts::NZ, tuple<omp_lock_t, vector<tuple<size_t, size_t>>>());
	for (size_t i = 0; i < consts::NX*consts::NZ; i++) {
		omp_init_lock(&get<0>(marked[i]));
	}
	return Beam {
		rays,
		marked,
		vector<tuple<bool, tuple<size_t, size_t>>>(),
	};
}

Beam beam2() {
	double dx = (consts::BEAM_MAX_Z-consts::BEAM_MIN_Z)/((double)consts::NRAYS-1.0);
	vector<Ray> rays;
	rays.reserve(consts::NRAYS);
	for (size_t i = 0; i < consts::NRAYS; i++) {
		double x0 = (double)i * dx + consts::BEAM_MIN_Z + consts::OFFSET2;
		rays.push_back(Ray {
			vector<Crossing>(), // crossings
			x0, // x0
			consts::ZMIN+0.1e-4, // z0
			x0 + consts::CHILD_OFFSET, // cx0
			consts::ZMIN+0.1e-4, // cz0
			consts::DIR2[0], // kx0
			consts::DIR2[1], // kz0
		});
	}
	vector<tuple<omp_lock_t, vector<tuple<size_t, size_t>>>> marked
		(consts::NX*consts::NZ, tuple<omp_lock_t, vector<tuple<size_t, size_t>>>());
	for (size_t i = 0; i < consts::NX*consts::NZ; i++) {
		omp_init_lock(&get<0>(marked[i]));
	}
	return Beam {
		rays,
		marked,
		vector<tuple<bool, tuple<size_t, size_t>>>(),
	};
}
