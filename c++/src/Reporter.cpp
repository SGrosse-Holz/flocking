/* This file is distributed under MIT license as part of the project flocksims.
 * See LICENSE file for details.
 * (c) Simon Grosse-Holz, 2019
 */

# include "Reporter.h"
# include "macros.h"

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

HDF5Reporter::savingConformation::savingConformation(const State& state, reportMode toReport)
	: Conformation(state.conf)
{
	t = state.t;
	this->toReport = toReport;

	// Never report thm, it's useless
	delete[] thm;
	thm = 0;

	if (toReport & reportModes::noParticles)
	{
		delete[] x;
		delete[] y;
		delete[] thp;
		
		x = 0;
		y = 0;
		thp = 0;
	}
	if (!(toReport & reportModes::noAnalysis))
	{
		P = state.calc_polarization();
		dSdt_discrete = state.calc_dSdt_discrete();
	}
	if (!(toReport & reportModes::noTheory))
	{
		dSdt_expected = state.calc_dSdt_expected();
	}
}

void HDF5Reporter::report(const State& state, reportMode toReport)
{
	if (!do_report) return;

	datalist.emplace_back(state, toReport);
}

// The workhorse of the reporter. Aggregates all the data in the data list into
// big arrays and then writes those to one HDF5 file. In doing so, creates a
// new group in the file, such that this function can be called whenever we
// want to dump stuff from memory.
void HDF5Reporter::dump()
{
	if (!do_report) return;

	// Pool all the datasets from the list in big arrays
	// Note that the particle and analysis variables have their own times
	// associated with them, because they might be different.
	// Note: bookkeeping is easier if we aggregate the particle data, i.e.
	//       x, y, thm, thp in one 3d array. It will still be
	//       written as 3 named 1d arrays.
	double **particle_data=0, *tps=0;
	std::string pd_column_names[] = {"x", "y", "theta+"};
	double *tas=0, *Ses=0, *Ps=0;
	double *tts=0, *Sts=0;

	// Get the relevant sizes
	int Tp=0, Ta=0, Tt=0;
	for (const auto& conf : datalist)
	{
		if (!(conf.toReport & reportModes::noParticles))
			Tp++;
		if (!(conf.toReport & reportModes::noAnalysis))
			Ta++;
		if (!(conf.toReport & reportModes::noTheory))
			Tt++;
	}

	// Allocate memory
	// Note: make sure that the allocated memory is contiguous, otherwise
	//       HDF5 will not work
	int N = datalist.front().N;
	particle_data = new double*[4];
	for (int i = 0; i < 3; i++)
		particle_data[i] = new double[Tp*N];
	tps = new double[Tp];
	
	tas = new double[Ta];
	Ses = new double[Ta];
	Ps = new double[Ta];

	tts = new double[Tt];
	Sts = new double[Tt];

	int ip=0, ia=0, it=0;
	for (const auto& conf : datalist)
	{
		if (!(conf.toReport & reportModes::noParticles))
		{
			for (int in = 0; in < N; in++)
			{
				particle_data[0][ip*N + in] = conf.x[in];
				particle_data[1][ip*N + in] = conf.y[in];
				particle_data[2][ip*N + in] = conf.thp[in];
			}
			tps[ip] = conf.t;
			ip++;
		}

		if (!(conf.toReport & reportModes::noAnalysis))
		{
			Ses[ia] = conf.dSdt_discrete;
			Ps[ia] = conf.P;
			tas[ia] = conf.t;
			ia++;
		}

		if (!(conf.toReport & reportModes::noTheory))
		{
			Sts[it] = conf.dSdt_expected;
			tts[it] = conf.t;
			it++;
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

	for (int i = 0; i < 3; i++)
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
	dataset = outfile.createDataSet(group_name+"/dSdt_discrete", datatype, dataspace);
	dataset.write(Ses, datatype);
	dataset = outfile.createDataSet(group_name+"/P+", datatype, dataspace);
	dataset.write(Ps, datatype);

	// And the theory stuff
	dims[0] = (hsize_t)Tt;
	dataspace = DataSpace(1, dims);
	dataset = outfile.createDataSet(group_name+"/t_theory", datatype, dataspace);
	dataset.write(tas, datatype);
	dataset = outfile.createDataSet(group_name+"/dSdt_expected", datatype, dataspace);
	dataset.write(Sts, datatype);

	// Clean up
	outfile.close();
	delete group;

	for (int i = 0; i < 3; i++)
		delete[] particle_data[i];
	delete[] particle_data;
	delete[] tps;

	delete[] tas;
	delete[] Ses;
	delete[] Ps;

	delete[] tts;
	delete[] Sts;
}
