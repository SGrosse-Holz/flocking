# include "Reporter.h"

# include <fstream>
# include <iostream>
# include <sstream>

using namespace flocksims;
using namespace H5;

BaseReporter::BaseReporter(std::string filename) : filename(filename), do_report(true)
{
	if (filename.empty())
	{
		do_report = false;
		return;
	}

	try
	{
		outfile = H5File(filename, H5F_ACC_TRUNC);
		outfile.close();
	}
	catch (FileIException)
	{
		std::cerr << "Could not open file " << filename << " for writing!" << std::endl;
		do_report = false;
	}
}

void HDF5Reporter::report(const State& state, reportMode toReport)
{
	if (!do_report) return;

	datalist.push_back(state);

	if (toReport & reportModes::noParticles)
	{
		delete[] datalist.back().xm;
		delete[] datalist.back().ym;
		delete[] datalist.back().xp;
		delete[] datalist.back().yp;
		delete[] datalist.back().theta;

		datalist.back().xm = 0;
		datalist.back().ym = 0;
		datalist.back().xp = 0;
		datalist.back().yp = 0;
		datalist.back().theta = 0;
	}
	if (toReport & reportModes::noAnalysis)
		datalist.back().polarization = -2;
	
}

// The workhorse of the reporter. Aggregates all the data in the data list into
// big arrays and then writes those to one HDF5 file. In doing so, creates a
// new group in the file, such that this function can be called whenever we
// want to dump stuff from memory.
void HDF5Reporter::dump()
{
	if (!do_report) return;

	// Pool all the datasets from the queue in big arrays
	// Note that the particle and analysis variables have their own times
	// associated with them, because they might be different.
	// Note: bookkeeping is easier, if we aggregate the particle data, i.e.
	//       xm, ym, xp, yp, theta in one 3d array. It will still be
	//       written as 5 named 1d arrays.
	double **particle_data=0, *tps=0;
	std::string pd_column_names[] = {"x-", "y-", "x+", "y+", "theta"};
	double *tas=0, *Ss=0, *Ps=0;

	// Get the relevant sizes
	int Tp=0, Ta=0;
	for (const auto& state : datalist)
	{
		if (state.xm && state.ym && state.xp && state.yp && state.theta)
			Tp++;
		if (state.polarization > -2)
			Ta++;
	}

	// Allocate memory
	// Note: make sure that the allocated memory is contiguous, otherwise
	//       HDF5 will not work
	int N = datalist.front().N;
	particle_data = new double*[5];
	for (int i = 0; i < 5; i++)
		particle_data[i] = new double[Tp*N];
	tps = new double[Tp];
	
	tas = new double[Ta];
	Ss = new double[Ta];
	Ps = new double[Ta];

	int ip=0, ia=0;
	for (const auto& state : datalist)
	{
		if (state.xm && state.ym && state.xp && state.yp && state.theta)
		{
			for (int in = 0; in < N; in++)
			{
				particle_data[0][ip*N + in] = state.xm[in];
				particle_data[1][ip*N + in] = state.ym[in];
				particle_data[2][ip*N + in] = state.xp[in];
				particle_data[3][ip*N + in] = state.yp[in];
				particle_data[4][ip*N + in] = state.theta[in];
			}
			tps[ip] = state.t;
			ip++;
		}

		if (state.polarization > -2)
		{
			Ss[ia] = state.dS_past;
			Ps[ia] = state.polarization;
			tas[ia] = state.t;
			ia++;
		}
	}

	// Now write these big arrays to an hdf5 file, creating a new group
	// every time we dump() (which should be roughly once per simulation)
	static int groupID = 0;

	std::stringstream ss;
	ss << "/" << groupID++;
	std::string group_name = ss.str();

	outfile.openFile(filename, H5F_ACC_RDWR);
	Group *group = new Group(outfile.createGroup(group_name));

	// Particle data
	hsize_t dims[] = {(hsize_t)Tp, (hsize_t)N};
	DataSpace dataspace(2, dims);
	FloatType datatype(PredType::NATIVE_DOUBLE);
	DataSet dataset;

	for (int i = 0; i < 5; i++)
	{
		dataset = outfile.createDataSet(group_name+"/"+pd_column_names[i], datatype, dataspace);
		dataset.write(particle_data[i], datatype);
	}
	// The corresponding times
	dataspace = DataSpace(1, dims); // dims[0] = Tp already
	dataset = outfile.createDataSet(group_name+"/t_particles", datatype, dataspace);
	dataset.write(tps, datatype);

	// Now the analysis data
	dims[0] = (hsize_t)Ta;
	dataspace = DataSpace(1, dims);
	dataset = outfile.createDataSet(group_name+"/t_analysis", datatype, dataspace);
	dataset.write(tas, datatype);
	dataset = outfile.createDataSet(group_name+"/dS", datatype, dataspace);
	dataset.write(Ss, datatype);
	dataset = outfile.createDataSet(group_name+"/P", datatype, dataspace);
	dataset.write(Ps, datatype);

	// Clean up
	outfile.close();
	delete group;

	for (int i = 0; i < 5; i++)
		delete[] particle_data[i];
	delete[] particle_data;
	delete[] tps;

	delete[] tas;
	delete[] Ss;
	delete[] Ps;
}
