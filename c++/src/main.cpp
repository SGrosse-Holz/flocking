/* This file is distributed under MIT license as part of the project flocksims.
 * See LICENSE file for details.
 * (c) Simon Grosse-Holz, 2019
 */

# include "flocksims.h"
# include "OptionHandler.h"

# include "H5Cpp.h"

# include <cmath>
# include <iostream>
# include <string>

using namespace flocksims;

int main(int argc, char *argv[])
{
	// Option handling
	op::OptionHandler OH(R"HELP(
This program runs a simulation vaguely resembling the Viscek model. This means
the equations of motion are given by

	x_new = x_old + v_0*dt*cos(theta)
	y_new = y_old + v_0*dt*sin(theta)

	theta_new = theta_old - dt*Sum_j J*exp(-r_ij/d)*sin(theta_old - theta_j) + sqrt(2*T*dt)*eta_i

where eta_i a normally distributed Gaussian variable, v_0, J, d, T are
parameters of the model. Note that units allow us to set v_0 = 1 and at
the same time (PBC boxsize) = 1 as well.

Note that the positions and orientations (theta) are evaluated in an
alternating (leapfrog) pattern. Thus the symmetric formulation of a
'system state' is x-, y-, theta, x+, y+, which is redundant (but symmetric).
)HELP");

	int N = 100;
	int analyzeEvery = 1;
	double dt = 0.001;
	double dt_save = 1;
	double total_time = 10.0;
	double T = 1.0;
	double J = 1.0;
	double d_int = 0.1;
	std::string outfile = "";

	OH.addOption((new op::SingleValueOption<int>("N", N))->description(
				"Total number of particles (birds)\n\tdefault: 100"));
	// OH.addOption((new op::SingleValueOption<int>("spb", stepsPerBlock))->description(
				// "Simulation steps per block\n\tdefault: 1000"));
	OH.addOption((new op::SingleValueOption<double>("dt_save", dt_save))->description(
				"Time step between saving full conformations.\n\t \
	Set to 0 to never save conformations.\n\t \
	default: 1"));
	OH.addOption((new op::SingleValueOption<int>("ae", analyzeEvery))->description(
				"AnalyzeEvery: how many integration steps\n\t \
	between logging of polarization and\n\t \
	entropy production.\n\t \
	default: 10"));
	// OH.addOption((new op::SingleValueOption<int>("n", total_blocks))->description(
				// "Total number of blocks\n\tdefault: 10"));
	OH.addOption((new op::SingleValueOption<double>("runtime", total_time))->description(
				"Total runtime of the simulation.\n\tdefault: 10"));
	OH.addOption((new op::SingleValueOption<double>("T", T))->description(
				"Temperature (noise amplitude)\n\tdefault: 1.0"));
	OH.addOption((new op::SingleValueOption<double>("J", J))->description(
				"Interaction strength\n\tdefault: 1.0"));
	OH.addOption((new op::SingleValueOption<double>("d", d_int))->description(
				"Interaction range (exponential falloff distance)\n\tdefault: 0.1"));
	OH.addOption((new op::SingleValueOption<double>("dt", dt))->description(
				"Integration timestep\n\tdefault: 0.001"));
	OH.addOption((new op::SingleValueOption<std::string>("o", outfile))->description(
				"Output file name (.h5)\n\tdefault: no output"));

	op::pRes res = OH.procOptions(argc, argv);
	if(res != op::ok) return res;
	// Post proc
	bool report_states = dt_save > 0;
	dt_save = (dt_save>0) ? dt_save : total_time;
	int stepsPerBlock = (int)ceil(dt_save/dt);
	int analysesPerBlock = ceil(stepsPerBlock/float(analyzeEvery));
	stepsPerBlock = analysesPerBlock*analyzeEvery;
	int total_blocks = (int)ceil(total_time/dt_save);

	// Output
	BaseReporter *reporter;
	reporter = new HDF5Reporter(outfile);

	reporter->report("N", N);
	reporter->report("stepsPerBlock", stepsPerBlock);
	reporter->report("T", T);
	reporter->report("J", J);
	reporter->report("interaction range", d_int);
	reporter->report("dt", dt);

	// Simulation
	State state(N);
	state.initialize_randomly();
	if (report_states)
		reporter->report(state, reportModes::noAnalysis);

	Integrator integrator(T, J, d_int, dt);

	System system(state, integrator);

	for(int i = 0; i < total_blocks; i++)
	{
		for (int j = 0; j < analysesPerBlock-1; j++)
		{
			system.step(analyzeEvery);
			reporter->report(system.get_state(), reportModes::noParticles);
		}
		system.step(analyzeEvery);
		reporter->report(system.get_state(), report_states ? reportModes::all : reportModes::noParticles);

		// system.print();
	}

	// Game over, clean up
	reporter->dump();
	delete reporter;
	return 0;
}
