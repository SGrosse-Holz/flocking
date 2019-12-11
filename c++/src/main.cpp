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
	int nsteps = 0;
	int nsteps_max = 0;
	bool doTheory = false;
	bool report_states = false;
	bool hbd = false;
	double v0 = 1.0;
	double density = 1;
	double dt = 0.001;
	double dt_save = 0; // Save every conformation (if saving at all)
	double total_time = 10.0;
	double T = 1.0;
	double J = 1.0;
	double d_int = 1;
	double adaptive_precision = 0.0;
	std::string outfile = "";

	OH.addOption((new op::SingleValueOption<int>("steps", nsteps))->description(
				"Number of steps to take. If specified, this overrides -runtime."));
	OH.addOption((new op::SingleValueOption<int>("max_steps", nsteps_max))->description(
				"Maximum number of steps. This is useful with adaptive_precision."));
	OH.addOption((new op::SingleValueOption<int>("N", N))->description(
				"Total number of particles (birds)\n\tdefault: 100"));
	OH.addOption((new op::SingleValueOption<double>("n", density))->description(
				"Density of birds\n\tdefault: 1"));
	OH.addOption((new op::SingleValueOption<double>("v0", v0))->description(
				"Velocity of the birds\n\tdefault: 1"));
	OH.addOption((new op::SingleValueOption<double>("dt_save", dt_save))->description(
				"Time step between saving full conformations.\n\t \
	Set to 0 to never save conformations.\n\t \
	default: 1"));
	OH.addOption((new op::SingleValueOption<int>("ae", analyzeEvery))->description(
				"AnalyzeEvery: how many integration steps\n\t \
	between logging of polarization and\n\t \
	entropy production.\n\t \
	default: 10"));
	OH.addOption((new op::FlagOption("doTheory", doTheory))->description(
				"If specified, calculate the expected entropy production for each analyzed state."));
	OH.addOption((new op::FlagOption("reportStates", report_states))->description(
				"If set, write some states out with the data. c.f. dt_save."));
	OH.addOption((new op::FlagOption("hbd", hbd))->description(
				"Start from hbd config, instead of randomly"));
	OH.addOption((new op::SingleValueOption<double>("runtime", total_time))->description(
				"Total runtime of the simulation.\n\tdefault: 10"));
	OH.addOption((new op::SingleValueOption<double>("T", T))->description(
				"Temperature (noise amplitude)\n\tdefault: 1.0"));
	OH.addOption((new op::SingleValueOption<double>("J", J))->description(
				"Interaction strength\n\tdefault: 1.0"));
	OH.addOption((new op::SingleValueOption<double>("d", d_int))->description(
				"Interaction range (exponential falloff distance)\n\tdefault: 1"));
	OH.addOption((new op::SingleValueOption<double>("dt", dt))->description(
				"Integration timestep\n\tdefault: 0.001"));
	OH.addOption((new op::SingleValueOption<double>("adaptive_precision", adaptive_precision))->description(
				"Precision of the adaptive time step.\n\t \
	Reasonable values are 1e-2, default 0, in which\n\t \
	case the timestep will not be adapted."));
	OH.addOption((new op::SingleValueOption<std::string>("o", outfile))->description(
				"Output file name (.h5)\n\tdefault: no output"));

	op::pRes res = OH.procOptions(argc, argv);
	if(res != op::ok) return res;

	// Adaptive timestep
	if (adaptive_precision > 0)
		dt = adaptive_precision*std::min(std::min( 1/(J*N),  adaptive_precision/(2*T) ),
						 std::min( d_int/v0, 1/(v0*sqrt(density))     ));

	// Post proc
	dt_save = (dt_save > dt) ? dt_save : dt;
	int stepsPerBlock = (int)ceil(dt_save/dt);
	int analysesPerBlock = ceil(stepsPerBlock/float(analyzeEvery));
	stepsPerBlock = analysesPerBlock*analyzeEvery;
	dt_save = dt*stepsPerBlock;
	int total_blocks = (int)ceil(total_time/dt_save);
	if (nsteps > 0)
		total_blocks = (int)ceil(nsteps/float(stepsPerBlock));
	if (nsteps_max > 0 && total_blocks*stepsPerBlock > nsteps_max)
		total_blocks = (int)ceil(nsteps_max/float(stepsPerBlock));
	total_time = total_blocks*dt_save;

	// Output
	BaseReporter *reporter;
	reporter = new HDF5Reporter(outfile);

	reporter->report("N", N);
	reporter->report("density", density);
	reporter->report("total runtime", total_time);
	reporter->report("T", T);
	reporter->report("J", J);
	reporter->report("interaction range", d_int);
	reporter->report("dt", dt);
	reporter->report("v0", v0);

	// Simulation
	DBGOUT("initializing...");
	Conformation start_conf(N, density);
	if (hbd)
	{
		try
		{
			start_conf.initialize_hbd();
		}
		catch (const char* msg)
		{
			std::cout << msg << std::endl;
			return 1;
		}
	}
	else
		start_conf.initialize_randomly();
	State initial_state(start_conf, T, J, dt, d_int, v0);
	DBGOUT("finished initializing");

	LeapFrogIntegrator integrator;
	System system(initial_state, &integrator);

	// system.get_state().dt = 1e-4;
	system.step(10000000);
	// system.get_state().dt = dt;

	if (report_states)
		reporter->report(system.get_state(), reportModes::noAnalysis | reportModes::noTheory);

	DBGOUT("starting run...");
	for(int i = 0; i < total_blocks; i++)
	{
		for (int j = 0; j < analysesPerBlock-1; j++)
		{
			DBGOUT("stepping...");
			system.step(analyzeEvery);
			system.report(reporter, reportModes::noParticles | (!doTheory * reportModes::noTheory));
		}
		DBGOUT("stepping...");
		system.step(analyzeEvery);
		system.report(reporter, (report_states ? reportModes::all : reportModes::noParticles)
				      | (!doTheory * reportModes::noTheory));

		// system.print();
	}

	// Game over, clean up
	reporter->dump();
	delete reporter;
	return 0;
}
