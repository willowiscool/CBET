#ifndef MESH_H
#define MESH_H

#include <vector>
#include <tuple>

using namespace std;

struct Point {
	double x;
	double z;
	double eden;
	double machnum;
	double kib_multiplier;
	double permittivity_multiplier;
};

struct Mesh {
	Point* points;
	double dx;
	double dz;
	size_t nx;
	size_t nz;
	double xmin;
	double xmax;
	double zmin;
	double zmax;

	static Mesh new_lin(double xmin, double xmax, size_t nx, double zmin, double zmax, size_t nz);

	Point& get(size_t x, size_t z);
	tuple<size_t, size_t> get_mesh_coords(double x, double z);
	tuple<size_t, size_t> get_mesh_coords_in_area(double x, double z, tuple<size_t, size_t> minpt, tuple<size_t, size_t> maxpt);
};

Mesh new_mesh();

#endif
