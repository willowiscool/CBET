#include <tuple>
#include <vector>

#include "cbet.hpp"
#include "mesh.hpp"
#include "beam.hpp"
#include "consts.hpp"

using namespace std;

void cbet(Mesh& mesh, vector<Beam>& beams) {
	double currmax = consts::MAX_INCR;
	create_raystore(beams, mesh.nx, mesh.nz);

	for (size_t i = 1; i <= 500; i++) {
		get_cbet_gain(mesh, beams);
		double updateconv = update_intensities(beams, 0.0, currmax);
		if (updateconv <= consts::CONVERGE) {
			break;
		}

		double currmaxa = consts::MAX_INCR*pow(consts::CBETCONVERGENCE, i);
		double currmaxb = consts::CBETCONVERGENCE*updateconv;
		currmax = min(currmaxa, currmaxb);
	}
	post(mesh, beams);
}

void post(Mesh& mesh, vector<Beam>& beams) {
	double w0 = 2.0*M_PI*consts::C_SPEED/consts::LAMBDA;
	double norm_factor_const = sqrt(8.0*M_PI/consts::C_SPEED) * consts::ESTAT / (consts::ME_G*consts::C_SPEED*w0) * sqrt(1e14*1e7);

	for (size_t beamnum = 0; beamnum < beams.size(); beamnum++) {
		for (size_t raynum = 0; raynum < beams[beamnum].rays.size(); raynum++) {
			Ray& ray = beams[beamnum].rays[raynum];
			for (size_t cnum = 0; cnum < ray.crossings.size(); cnum++) {
				Crossing& crossing = ray.crossings[cnum];
				Crossing& crossing_next = ray.crossings[min(cnum+1, ray.crossings.size()-1)];
				size_t ix = crossing.boxesx;
				size_t iz = crossing.boxesz;

				double area_avg = (crossing.area_ratio+crossing_next.area_ratio)/2.0;
				double ne_over_nc = mesh.get(ix, iz).eden/consts::NCRIT;
				double ne_over_nc_corrected = min(ne_over_nc, 1.0);
				double ne_term = sqrt(1.0-ne_over_nc_corrected);
				double epsilon_eff = ne_term*ne_term;
				double interaction_mult = 1.0/(area_avg*ne_term)*1.0/sqrt(epsilon_eff);
				double norm_factor = norm_factor_const * sqrt(interaction_mult) * pow(epsilon_eff, 0.25);
				double intensity_new = sqrt(crossing.i_b) * norm_factor;
				
				// now we have mutable borrow
				crossing.i_b = intensity_new;
			}
		}
	}
}

double update_intensities(vector<Beam>& beams, double conv_max, double curr_max) {
	double curr_conv_max = conv_max;
	for (size_t beamnum = 0; beamnum < beams.size(); beamnum++) {
#pragma omp parallel for reduction(max: curr_conv_max)
		for (size_t raynum = 0; raynum < beams[beamnum].rays.size(); raynum++) {
			Ray& ray = beams[beamnum].rays[raynum];
			double i0 = ray.crossings[0].i_b;
			double mult_acc = 1.0;
			for (size_t cnum = 0; cnum < ray.crossings.size(); cnum++) {
				Crossing& crossing = ray.crossings[cnum];
				auto [new_intensity, new_conv_max] = limit_energy(crossing, mult_acc, i0, curr_max, curr_conv_max);
				curr_conv_max = new_conv_max;
				mult_acc *= crossing.w_mult;
				crossing.i_b = new_intensity;
			}
		}
	}
	return curr_conv_max;
}

tuple<double, double> limit_energy(Crossing& crossing, double multiplier_acc, double i0, double curr_max, double max_change) {
	double i_prev = crossing.i_b;
	double i_curr = i0*multiplier_acc;
	double fractional_change = abs(i_curr-i_prev)/i_prev;
	double new_max_change = max(fractional_change, max_change);

	if (fractional_change > curr_max) {
		double sign = i_curr - i_prev > 0.0 ? 1.0 : -1.0;
		double correction = 1.0 + curr_max*sign;
		i_curr = i_prev*correction;
	}
	return {i_curr, new_max_change}; 
}

void get_cbet_gain(Mesh& mesh, vector<Beam>& beams) {
	for (size_t beamnum = 0; beamnum < beams.size(); beamnum++) {
#pragma omp parallel for
		for (size_t raynum = 0; raynum < beams[beamnum].rays.size(); raynum++) {
			for (size_t cnum = 0; cnum < beams[beamnum].rays[raynum].crossings.size(); cnum++) {
				Crossing& crossing = beams[beamnum].rays[raynum].crossings[cnum];
				size_t ix = crossing.boxesx;
				size_t iz = crossing.boxesz;

				auto other_beam_cbet_incr = [&] (Beam& other_beam) {
					// raystore value
					auto [has_crossing, inds] = other_beam.raystore[ix*mesh.nz+iz];
					// tuple<bool, tuple<size_t, size_t>> rval = other_beam.raystore[ix*mesh.nz+iz];
					if (!has_crossing) {
						return 0.0;
					}
					// giving up on not using auto (no idea why I haven't been? whatever)
					vector<Crossing>& other_ray_crossings = other_beam.rays[get<0>(inds)].crossings;
					Crossing& raycross = other_ray_crossings[get<1>(inds)];
					Crossing& raycross_next =
						get<1>(inds)+1 == other_ray_crossings.size()
						? raycross
						: other_ray_crossings[get<1>(inds)+1];

					return get_cbet_increment(mesh, crossing, raycross, raycross_next);
				};

				double cbet_sum = 0;
				for (size_t o_bnum = 0; o_bnum < beams.size(); o_bnum++) {
					if (o_bnum == beamnum) continue;
					cbet_sum += other_beam_cbet_incr(beams[o_bnum]);
				}
				crossing.w_mult = exp(-1.0*cbet_sum);
			}
		}
	}
}

double get_cbet_increment(Mesh& mesh, Crossing& crossing, Crossing& raycross, Crossing& raycross_next) {
	Point& mesh_pt = mesh.get(crossing.boxesx, crossing.boxesz);
	// comments been moved between implementations like thrice honestly lol
	
	// INTERACTION MULTIPLIER: find ray_o's interaction multiplier
	// get the avereage area of the ray across the zone
	double area_avg = (raycross.area_ratio+raycross_next.area_ratio)/2.0;
	// NOTE: This neOverNc value can be taken straight from the grid
	double ne_over_nc = mesh_pt.eden/consts::NCRIT;
	// clamp the maximum neOverNc value to
	double ne_over_nc_corrected = min(ne_over_nc, 1.0);
	// term used multiple times in calculations
	double ne_term = sqrt(1.0-ne_over_nc_corrected);
	double epsilon_eff = ne_term*ne_term;
	double interaction_mult = 1.0/(area_avg*ne_term)*1.0/sqrt(epsilon_eff);
	
	// eta Routine
	double kx_seed = crossing.dkx;
	double kz_seed = crossing.dkz;
	double kx_pump = raycross.dkx;
	double kz_pump = raycross.dkz;
	
	double machx = mesh_pt.machnum; // the mach number of the plasma velocity
	double machz = 0.0; // TODO: Implement multidimensional plasma velocity
	
	// assuming uniform wavelength/frequency
	double omega1 = consts::OMEGA;
	double omega2 = consts::OMEGA;
	// find ion acoustic wave vector, difference of ray trajectories scaled to
	// omega/(sqrt(1-neOverNc)*c)
	
	double iaw_vector[] = {
	    (omega1*kx_seed - omega2*kx_pump)*sqrt(1.0-ne_over_nc)/consts::C_SPEED,
	    (omega1*kz_seed - omega2*kz_pump)*sqrt(1.0-ne_over_nc)/consts::C_SPEED,
	};
	double k_iaw = sqrt(iaw_vector[0]*iaw_vector[0] + iaw_vector[1]*iaw_vector[1]); // magnitude of iaw vector
	double eta_numerator = omega1-omega2 - (iaw_vector[0]*machx + iaw_vector[1]*machz)*consts::CS;
	double eta_denominator = k_iaw*consts::CS;
	double eta = eta_numerator/eta_denominator;
	
	// FIND COUPLING MULTIPLIER
	// Split coupling multiplier into discrete chuncks [sic]->easier debugging
	double param1 = consts::CBET_CONST/(consts::OMEGA*(consts::TE_EV/1e3 + 3.0 * consts::TI_EV/1e3/consts::Z));
	double param2 = ne_over_nc/consts::IAW*consts::IAW*consts::IAW*eta; // Need to fix iaw
	double param3 = pow(eta*eta-1.0, 2) + consts::IAW*consts::IAW*eta*eta;
	double param4 = interaction_mult;
	// WILLOW SAYS: "ds" is replaced with crossing.dkmag, cuz that's what it is
	double coupling_mult = param1*param2/param3*param4*crossing.dkmag; // get the coupling multiplier
	
	double other_intensity1 = raycross.i_b;
	double other_intensity2 = raycross_next.i_b;
	// average the intensity of the other ray across the cell
	double avg_intensity = (other_intensity1+other_intensity2)/2.0;
	
	return coupling_mult*avg_intensity;
}

tuple<Crossing&, Crossing&> get_raycross(Ray& ray, size_t ix, size_t iz) {
	size_t raycross = 0;
	while (ray.crossings[raycross].boxesx != ix || ray.crossings[raycross].boxesz != iz) {
		raycross += 1;
		// erasing panic (she is taking risks for an amazing reason)
	}
	return {ray.crossings[raycross], ray.crossings[min(ray.crossings.size()-1, raycross+1)]};
}

void init_crossings(vector<Beam>& beams, double intensity) {
	for (size_t beamnum = 0; beamnum < beams.size(); beamnum++) {
		double increment = (consts::BEAM_MAX_Z-consts::BEAM_MIN_Z)/((double)beams[beamnum].rays.size() - 1.0);
		for (size_t raynum = 0; raynum < beams[beamnum].rays.size(); raynum++) {
			double offset = consts::BEAM_MIN_Z + (increment*(double)raynum);
			double iintensity = (intensity/1e14)*exp(-2.0*pow(abs(offset/2e-4),4));
			// okay this is a lot of indirection haha girl just get a reference
			// whateverrrr maybe compiler optimizes it?
			for (size_t cnum = 0; cnum < beams[beamnum].rays[raynum].crossings.size(); cnum++) {
				beams[beamnum].rays[raynum].crossings[cnum].i_b = iintensity;
				beams[beamnum].rays[raynum].crossings[cnum].w_mult = 1.0;
			}
		}
	}
}

void create_raystore(vector<Beam>& beams, size_t nx, size_t nz) {
	for (size_t beamnum = 0; beamnum < beams.size(); beamnum++) {
		for (size_t i = 0; i < nx*nz; i++) {
			size_t x = i / nz;
			size_t z = i % nz;
			// TODO locks??
			vector<tuple<size_t, size_t>>& marked = get<1>(beams[beamnum].marked[x*nx+z]);
			if (marked.size() == 0) {
				beams[beamnum].raystore.push_back({false, {0, 0}});
			} else if (marked.size() == 1) {
				beams[beamnum].raystore.push_back({true, marked[0]});
			} else {
				beams[beamnum].raystore.push_back({true, marked[rand()%marked.size()]});
			}
		}
	}
}
