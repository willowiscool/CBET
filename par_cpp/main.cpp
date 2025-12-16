#include <iostream>
#include <vector>
#include <H5Cpp.h>
#include "ray_trace.hpp"
#include "beam.hpp"
#include "cbet.hpp"
#include "consts.hpp"

void save_hdf5(Mesh& mesh, vector<Beam>& beams, std::string filename);

int main() {
	cout << "Creating and initializing mesh" << endl;
	Mesh m = new_mesh();
	cout << "Creating and initializing beams" << endl;
	Beam b1 = beam1();
	Beam b2 = beam2();
	vector<Beam> beams = {b1, b2};
	cout << "Tracing rays" << endl;
	ray_trace(m, beams);
	cout << "Doing CBET calculation" << endl;
	init_crossings(beams, consts::INTENSITY);
	cbet(m, beams);
	cout << "Saving to \"out.hdf5\"" << endl;
	save_hdf5(m, beams, "out.hdf5");
	return 0;
}

void save_hdf5(Mesh& mesh, vector<Beam>& beams, std::string filename) {
	double* wplot = new double[mesh.nx*mesh.nz]();

	for (size_t beamnum = 0; beamnum < beams.size(); beamnum++) {
		// makes sure they're zeroed
		// but also like this is pretty not efficient?
		// but it happens only once per beam--it's fine
		double* intensities = new double[mesh.nx*mesh.nz]();
		size_t* ct = new size_t[mesh.nx*mesh.nz]();
		for (size_t rnum = 0; rnum < beams[beamnum].rays.size(); rnum++) {
			for (size_t cnum = 0; cnum < beams[beamnum].rays[rnum].crossings.size(); cnum++) {
				Crossing& c = beams[beamnum].rays[rnum].crossings[cnum];
				// cout << beamnum << "\t" << c.boxesx << "\t" << c.boxesz << endl;
				intensities[c.boxesx*mesh.nz + c.boxesz] += c.i_b;
				ct[c.boxesx*mesh.nz + c.boxesz] += 1;
			}
		}
		for (size_t x = 0; x < mesh.nx; x++) {
			for (size_t z = 0; z < mesh.nz; z++) {
				wplot[x*mesh.nz + z] += pow(intensities[x*mesh.nz + z] / max(ct[x*mesh.nz + z], (size_t)1), 2);
			}
		}
		delete[] intensities;
		delete[] ct;
	}

	for (size_t x = 0; x < mesh.nx; x++) {
		for (size_t z = 0; z < mesh.nz; z++) {
			wplot[x*mesh.nz + z] = sqrt(wplot[x*mesh.nz + z]);
		}
	}

	H5::H5File file(filename, H5F_ACC_TRUNC);
	H5::IntType datatype(H5::PredType::NATIVE_DOUBLE);
	datatype.setOrder(H5T_ORDER_LE);
	hsize_t dims[2] = {mesh.nx, mesh.nz};
	H5::DataSpace dataspace(2, dims);
	H5::DataSet dataset = file.createDataSet("/wplot", datatype, dataspace);
	dataset.write(wplot, H5::PredType::NATIVE_DOUBLE);

	delete[] wplot;
}
