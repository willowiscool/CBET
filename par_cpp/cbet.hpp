#ifndef CBET_H
#define CBET_H

#include <tuple>
#include <vector>
#include "mesh.hpp"
#include "beam.hpp"

using namespace std;

void cbet(Mesh& mesh, vector<Beam>& beams);
void post(Mesh& mesh, vector<Beam>& beams);
double update_intensities(vector<Beam>& beams, double conv_max, double curr_max);
tuple<double, double> limit_energy(Crossing& crossing, double multiplier_acc, double i0, double curr_max, double max_change);
void get_cbet_gain(Mesh& mesh, vector<Beam>& beams);
double get_cbet_increment(Mesh& mesh, Crossing& crossing, Crossing& raycross, Crossing& raycross_next);
tuple<Crossing&, Crossing&> get_raycross(Ray& ray, size_t ix, size_t iz);
void init_crossings(vector<Beam>& beams, double intensity);
void create_raystore(vector<Beam>& beams, size_t nx, size_t nz);

#endif
