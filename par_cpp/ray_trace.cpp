#include <cstddef>
#include <cassert>
#include <iostream>
#include <vector>
#include <tuple>
#include <omp.h>
#include "consts.hpp"
#include "beam.hpp"
#include "mesh.hpp"
#include "ray_trace.hpp"

using namespace std;

void ray_trace(Mesh& mesh, vector<Beam>& beams) {
	vector<tuple<double, double>> deden(mesh.nx*mesh.nz);
	for (size_t x = 0; x < mesh.nx - 1; x++) {
		for (size_t z = 0; z < mesh.nz - 1; z++) {
			Point& p = mesh.get(x, z);
			double dedendx = x == 0 ? 0 :
				(p.eden - mesh.get(x-1, z).eden)/(p.x-mesh.get(x-1, z).x);
			double dedendz = z == 0 ? 0 :
				(p.eden - mesh.get(x, z-1).eden)/(p.z-mesh.get(x, z-1).z);
			deden[x*mesh.nz+z] = {dedendx, dedendz};
		}
	}

	for (size_t i = 0; i < mesh.nz; i++) {
		std::get<0>(deden[i]) = std::get<0>(deden[mesh.nz+i]);
		deden[(mesh.nx-1)*mesh.nz+i] = deden[(mesh.nx-2)*mesh.nz+i];
	}
	for (size_t i = 0; i < mesh.nx; i++) {
		std::get<1>(deden[i*mesh.nz]) = std::get<1>(deden[i*mesh.nz+1]);
		deden[i*mesh.nz + (mesh.nz-1)] = deden[i*mesh.nz + (mesh.nz-1)];
	}

	for (size_t beamnum = 0; beamnum < beams.size(); beamnum++) {
#pragma omp parallel for
		for (size_t raynum = 0; raynum < beams[beamnum].rays.size(); raynum++) {
			// I know making and destroying a new vector is MASSIVELY inefficient
			// but it matches in both implementations so...
			tuple<vector<double>, vector<double>> child =
				launch_child_ray(beams[beamnum].rays[raynum], mesh, deden);
			launch_parent_ray(beams[beamnum].rays[raynum], mesh, deden,
				child, beams[beamnum].marked, raynum);
		}
	}
}

void launch_parent_ray(
	Ray& ray, Mesh& mesh, vector<tuple<double, double>>& deden,
	tuple<vector<double>, vector<double>> child,
	vector<tuple<omp_lock_t, vector<tuple<size_t, size_t>>>>& marked,
	size_t raynum
) {
	vector<double> childx = std::get<0>(child);
	vector<double> childz = std::get<1>(child);

	double x = ray.x0;
	double z = ray.z0;

	vector<double> distance;
	distance.reserve(childx.size());
	distance.push_back(0.0);
	for (size_t i = 1; i < childx.size(); i++) {
		distance.push_back(distance[i-1] + sqrt(pow(childx[i] - childx[i-1], 2) + pow(childz[i] - childz[i-1], 2)));
	}

	double init_diff_x = ray.x0 - childx[0];
	double init_diff_z = ray.z0 - childz[0];
	double init_diff_mag = sqrt(init_diff_x*init_diff_x + init_diff_z*init_diff_z);
	double init_proj_coeff = abs(ray.kx0 * (init_diff_z/init_diff_mag) - ray.kz0 * (init_diff_x/init_diff_mag));
	double init_area = init_diff_mag*init_proj_coeff;

	auto [meshx, meshz] = mesh.get_mesh_coords(ray.x0, ray.z0);
	double k = get_k(mesh, ray.x0);
	double knorm = sqrt(ray.kx0*ray.kx0 + ray.kz0*ray.kz0);
	double vx = consts::C_SPEED*consts::C_SPEED * ((ray.kx0 / knorm) * k) / consts::OMEGA;
	double vz = consts::C_SPEED*consts::C_SPEED * ((ray.kz0 / knorm) * k) / consts::OMEGA;

	double curr_dist = 0.0;
	for (size_t _ = 1; _ < consts::NT; _++) {
		double mydedendx = interp_closure(
			[&] (size_t i) -> double { return std::get<0>(deden[mesh.nz*i + (mesh.nz+1)/2]); },
			[&] (size_t i) -> double { return mesh.points[mesh.nz*i].x; }, mesh.nx,
			x
		);

		double prev_vx = vx;
		double prev_vz = vz;
		vx = vx - consts::C_SPEED*consts::C_SPEED / (2.0*consts::NCRIT) * mydedendx * consts::DT;
		vz = vz - consts::C_SPEED*consts::C_SPEED / (2.0*consts::NCRIT) * std::get<1>(deden[meshx*mesh.nz+meshz]) * consts::DT;
		double prev_x = x;
		double prev_z = z;
		x = prev_x + vx * consts::DT;
		z = prev_z + vz * consts::DT;

		size_t prev_meshx = meshx;
		size_t prev_meshz = meshz;
		tuple<size_t, size_t> mcoords = mesh.get_mesh_coords_in_area(
			x, z,
			{ meshx == 0 ? 0 : meshx-1, meshz == 0 ? 0 : meshz-1 },
			{ min(mesh.nx-1, meshx+1), min(mesh.nz-1, meshz+1) }
		);
		meshx = get<0>(mcoords);
		meshz = get<1>(mcoords);

		double lastx = 10000.0;
		double lastz = 10000.0;
		bool is_cross_x = false;
		bool is_cross_z = false;

		if (meshx != prev_meshx) {
			double currx = mesh.get(meshx, 0).x;
			if (!((x > currx && prev_x <= currx) || (x < currx && prev_x >= currx))) {
				currx = mesh.get(prev_meshx, 0).x;
			}
			// double crossx = interp({prev_z, z}, {prev_x, x}, currx);
			// double crossx = currx <= prev_x ? prev_z :
			// 	currx >= x ? z :
			// 	prev_z + ((z - prev_z) / (x - prev_x) * (currx - prev_x));
			 double crossx = prev_z + ((z - prev_z) / (x - prev_x) * (currx - prev_x)); // ***changed***
			double frac = (currx - prev_x) / (x - prev_x);

			if (abs(crossx - lastz) > 1.0e-12) {
				double extra = sqrt(pow(currx - prev_x, 2) + pow(crossx - prev_z, 2));
				double childxp = interp(childx, distance, curr_dist + extra);
				double childzp = interp(childz, distance, curr_dist + extra);

				double diff_x = currx - childxp;
				double diff_z = crossx - childzp;
				double diff_mag = sqrt(diff_x*diff_x + diff_z*diff_z);

				double interpkx = frac*prev_vx + (1.0-frac)*vx;
				double interpkz = frac*prev_vz + (1.0-frac)*vz;
				double interpk_mag = sqrt(interpkx*interpkx + interpkz*interpkz);

				double proj_coeff = abs((interpkx/interpk_mag)*(diff_z/diff_mag) - (interpkz/interpk_mag)*(diff_x/diff_mag));
				// Crossing
				ray.crossings.push_back({
					currx, // x
					crossx, // z
					meshx, // boxesx
					meshz, // boxesz
					diff_mag*proj_coeff/init_area, // area_ratio
					0.0, // dkx
					0.0, // dkz
					0.0, // dkmag
					-1.0, // i_b
					1.0, // w_mult
				});
				omp_set_lock(&get<0>(marked[meshx*mesh.nx + meshz]));
				get<1>(marked[meshx*mesh.nx + meshz]).push_back({raynum, ray.crossings.size()-1});
				omp_unset_lock(&get<0>(marked[meshx*mesh.nx + meshz]));

				is_cross_x = true;
				lastx = currx;
			}
		}
		if (meshz != prev_meshz) {
			double currz = mesh.get(0, meshz).z;
			if (!((z > currz && prev_z <= currz) || (z < currz && prev_z >= currz))) {
				currz = mesh.get(0, prev_meshz).z;
			}

			double crossz = prev_x + ((x - prev_x) / (z - prev_z)) * (currz - prev_z);
			double frac = (currz - prev_z) / (z - prev_z);

			if (abs(crossz - lastx) > 1.0e-12) {
				double extra = sqrt(pow(crossz - prev_x, 2) + pow(currz - prev_z, 2));
				double childxp = interp(childx, distance, curr_dist + extra);
				double childzp = interp(childz, distance, curr_dist + extra);

				double diff_x = crossz - childxp;
				double diff_z = currz - childzp;
				double diff_mag = sqrt(diff_x*diff_x + diff_z*diff_z);

				double interpkx = frac*prev_vx + (1.0-frac)*vx;
				double interpkz = frac*prev_vz + (1.0-frac)*vz;
				double interpk_mag = sqrt(interpkx*interpkx + interpkz*interpkz);

				double proj_coeff = abs((interpkx/interpk_mag)*(diff_z/diff_mag) - (interpkz/interpk_mag)*(diff_x/diff_mag));
				// Crossing
				ray.crossings.push_back({
					crossz, // x
					currz, // z
					meshx, // boxesx
					meshz, // boxesz
					diff_mag*proj_coeff/init_area, // area_ratio
					0.0, // dkx
					0.0, // dkz
					0.0, // dkmag
					-1.0, // i_b
					1.0, // w_mult
				});
				omp_set_lock(&get<0>(marked[meshx*mesh.nx + meshz]));
				get<1>(marked[meshx*mesh.nx + meshz]).push_back({raynum, ray.crossings.size()-1});
				omp_unset_lock(&get<0>(marked[meshx*mesh.nx + meshz]));

				is_cross_z = true;
				lastz = currz;
			}
		}
		if (is_cross_x && is_cross_z) {
			size_t last_crossing_ind = ray.crossings.size()-1;
			if ((x - prev_x) * (ray.crossings[last_crossing_ind].x - ray.crossings[last_crossing_ind-1].x) < 0.0) {
				swap(ray.crossings[last_crossing_ind], ray.crossings[last_crossing_ind-1]);
			}
			if (last_crossing_ind > 1)
				calc_dk(ray, last_crossing_ind-2);
			calc_dk(ray, last_crossing_ind-1);
		} else if ((is_cross_x || is_cross_z) && ray.crossings.size() > 1) {
			calc_dk(ray, ray.crossings.size()-2);
		}

		curr_dist += sqrt(pow(x - prev_x, 2) + pow(z - prev_z, 2));
		if (x < mesh.xmin || x > mesh.xmax || z < mesh.zmin || z > mesh.zmax) {
			break;
		}
	}
}
void calc_dk(Ray& ray, size_t ind) {
	double dkx = ray.crossings[ind+1].x - ray.crossings[ind].x;
	double dkz = ray.crossings[ind+1].z - ray.crossings[ind].z;
	double dkmag = sqrt(pow(dkx, 2) + pow(dkz, 2));
	double dkx_new = dkx/dkmag;
	double dkz_new = dkz/dkmag;
	double dkmag_new = dkmag*10000.0;

	ray.crossings[ind].dkx = dkx_new;
	ray.crossings[ind].dkz = dkz_new;
	ray.crossings[ind].dkmag = dkmag_new;
}

tuple<vector<double>, vector<double>> launch_child_ray(
	Ray& ray, Mesh& mesh, vector<tuple<double, double>>& deden
) {
	vector<double> x;
	vector<double> z;

	x.push_back(ray.cx0);
	z.push_back(ray.cz0);

	auto [meshx, meshz] = mesh.get_mesh_coords(ray.cx0, ray.cz0);
	double k = get_k(mesh, ray.cx0);
	double knorm = sqrt(ray.kx0*ray.kx0 + ray.kz0*ray.kz0);
	double vx = consts::C_SPEED*consts::C_SPEED * ((ray.kx0 / knorm) * k) / consts::OMEGA;
	double vz = consts::C_SPEED*consts::C_SPEED * ((ray.kz0 / knorm) * k) / consts::OMEGA;

	for (size_t tt = 1; tt < consts::NT; tt++) {
		// calc. vel & pos
		// assumes linear gradient w/ constant mydedendz
		double mydedendx = interp_closure(
			[&] (size_t i) -> double { return std::get<0>(deden[mesh.nz*i + (mesh.nz+1)/2]); },
			[&] (size_t i) -> double { return mesh.points[mesh.nz*i].x; }, mesh.nx,
			x[tt-1]
		);

		vx = vx - consts::C_SPEED*consts::C_SPEED / (2.0*consts::NCRIT) * mydedendx * consts::DT;
		vz = vz - consts::C_SPEED*consts::C_SPEED / (2.0*consts::NCRIT) * std::get<1>(deden[meshx*mesh.nz+meshz]) * consts::DT;

		x.push_back(x[tt-1]+vx * consts::DT);
		z.push_back(z[tt-1]+vz * consts::DT);

		tuple<double, double> mcoords = mesh.get_mesh_coords_in_area(
			x[tt], z[tt],
			{ meshx == 0 ? 0 : meshx-1, meshz == 0 ? 0 : meshz-1 },
			{ min(mesh.nx-1, meshx+1), min(mesh.nz-1, meshz+1) }
		);
		meshx = std::get<0>(mcoords);
		meshz = std::get<1>(mcoords);

		if (x[tt] < mesh.xmin || x[tt] > mesh.xmax || z[tt] < mesh.zmin || z[tt] > mesh.zmax) break;
	}

	return {x, z};
}

double get_k(Mesh& mesh, double x0) {
	double wpe_interp = sqrt(
		interp_closure(
			[&] (size_t i) -> double { return mesh.points[i*mesh.nz].eden; },
			[&] (size_t i) -> double { return mesh.points[i*mesh.nz].x; }, mesh.nx,
			x0
		) * 1e6 * pow(consts::EC, 2) / (consts::ME*consts::E0)
	);

	return sqrt((pow(consts::OMEGA, 2) - pow(wpe_interp, 2)) / (pow(consts::C_SPEED, 2)));
}
// implementation of interp that takes closures, for get_k
template<typename Func1, typename Func2>
double interp_closure(Func1 y, Func2 x, size_t len, double xp) {
	if (x(0) <= x(len-1)) {
		// x monotonically increase
		if (xp <= x(0)) {
			return y(0);
		} else if (xp >= x(len-1)) {
			return y(len-1);
		}

		size_t low = 0;
		size_t high = len - 1;
		size_t mid = (low + high) >> 1;
		while (low < high - 1) {
			if (x(mid) >= xp) {
				high = mid;
			} else {
				low = mid;
			}
			mid = (low + high) >> 1;
		}

		assert((xp >= x(mid)) && (xp <= x(mid+1)));
		return y(mid) + ((y(mid+1) - y(mid)) / (x(mid+1) - x(mid)) * (xp - x(mid)));
	} else {
		if (xp >= x(0)) {
			return y(0);
		} else if (xp <= x(len-1)) {
			return y(len-1);
		}

		size_t low = 0;
		size_t high = len - 1;
		size_t mid = (low + high) >> 1;
		while (low < high - 1) {
			if (x(mid) <= xp) {
				low = mid;
			} else {
				high = mid;
			}
			mid = (low + high) >> 1;
		}

		assert((xp <= x(mid)) && (xp >= x(mid+1)));
		return y(mid) + ((y(mid+1)-y(mid))/(x(mid+1)-x(mid)) * (xp-x(mid)));
	}
}
double interp(vector<double>& y, vector<double>& x, double xp) {
	assert(y.size() == x.size());
	size_t len = x.size();
	if (x[0] <= x[len-1]) {
		// x monotonically increase
		if (xp <= x[0]) {
			return y[0];
		} else if (xp >= x[len-1]) {
			return y[len-1];
		}

		size_t low = 0;
		size_t high = len - 1;
		size_t mid = (low + high) >> 1;
		while (low < high - 1) {
			if (x[mid] >= xp) {
				high = mid;
			} else {
				low = mid;
			}
			mid = (low + high) >> 1;
		}

		assert((xp >= x[mid]) && (xp <= x[mid+1]));
		return y[mid] + ((y[mid+1] - y[mid]) / (x[mid+1] - x[mid]) * (xp - x[mid]));
	} else {
		if (xp >= x[0]) {
			return y[0];
		} else if (xp <= x[len-1]) {
			return y[len-1];
		}

		size_t low = 0;
		size_t high = len - 1;
		size_t mid = (low + high) >> 1;
		while (low < high - 1) {
			if (x[mid] <= xp) {
				low = mid;
			} else {
				high = mid;
			}
			mid = (low + high) >> 1;
		}

		assert((xp <= x[mid]) && (xp >= x[mid+1]));
		return y[mid] + ((y[mid+1]-y[mid])/(x[mid+1]-x[mid]) * (xp-x[mid]));
	}
}
