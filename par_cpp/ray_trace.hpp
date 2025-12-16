#ifndef RAY_TRACE_H
#define RAY_TRACE_H

#include <vector>
#include <omp.h>
#include "mesh.hpp"
#include "beam.hpp"

void ray_trace(Mesh& mesh, vector<Beam>& beams);
tuple<vector<double>, vector<double>> launch_child_ray(
	Ray& ray, Mesh& mesh, vector<tuple<double, double>>& deden
);
void launch_parent_ray(
	Ray& ray, Mesh& mesh, vector<tuple<double, double>>& deden,
	tuple<vector<double>, vector<double>> child,
	vector<tuple<omp_lock_t, vector<tuple<size_t, size_t>>>>& marked,
	size_t raynum
);
double get_k(Mesh& mesh, double x0);
template<typename Func1, typename Func2>
double interp_closure(Func1 y, Func2 x, size_t len, double xp);
double interp(vector<double>& y, vector<double>& x, double xp);
void calc_dk(Ray& ray, size_t ind);

#endif
