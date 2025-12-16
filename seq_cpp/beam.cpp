#include <cstddef>
#include "beam.hpp"
#include "consts.hpp"

Beam beam1() {
	double dz = (consts::BEAM_MAX_Z-consts::BEAM_MIN_Z)/((double)consts::NRAYS-1.0);
	std::vector<Ray> rays;
	rays.reserve(consts::NRAYS);
	for (size_t i = 0; i < consts::NRAYS; i++) {
		double z0 = (double)i * dz + consts::BEAM_MIN_Z + consts::OFFSET1;
		rays.push_back(Ray {
			std::vector<Crossing>(), // crossings
			consts::XMIN, // x0
			z0,
			consts::XMIN, // cx0
			z0 + consts::CHILD_OFFSET, // cz0
			consts::DIR1[0], // kx0
			consts::DIR1[1], // kz0
		});
	}
	return Beam {
		rays,
		std::vector<std::vector<std::tuple<size_t, size_t>>>(consts::NX*consts::NZ, std::vector<std::tuple<size_t, size_t>>()),
		std::vector<std::tuple<bool, std::tuple<size_t, size_t>>>(),
	};
}

Beam beam2() {
	double dx = (consts::BEAM_MAX_Z-consts::BEAM_MIN_Z)/((double)consts::NRAYS-1.0);
	std::vector<Ray> rays;
	rays.reserve(consts::NRAYS);
	for (size_t i = 0; i < consts::NRAYS; i++) {
		double x0 = (double)i * dx + consts::BEAM_MIN_Z + consts::OFFSET2;
		rays.push_back(Ray {
			std::vector<Crossing>(), // crossings
			x0, // x0
			consts::ZMIN+0.1e-4, // z0
			x0 + consts::CHILD_OFFSET, // cx0
			consts::ZMIN+0.1e-4, // cz0
			consts::DIR2[0], // kx0
			consts::DIR2[1], // kz0
		});
	}
	return Beam {
		rays,
		std::vector<std::vector<std::tuple<size_t, size_t>>>(consts::NX*consts::NZ, std::vector<std::tuple<size_t, size_t>>()),
		std::vector<std::tuple<bool, std::tuple<size_t, size_t>>>(),
	};
}
