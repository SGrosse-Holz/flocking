/* This file is distributed under MIT license as part of the project flocksims.
 * See LICENSE file for details.
 * (c) Simon Grosse-Holz, 2019
 */

# include "Integrator.h"
# include "mathutils.h"
# include "macros.h"

# include <chrono>
# include <cmath>
# include <functional>
# include <random>
# include <iostream>

using namespace flocksims;

void LeapFrogIntegrator::step_orientations(State& state)
{
	// Status when this function is called:
	//  - state.conf.thm are the orientations we know
	//  - state.conf.x,y are already updated
	//  - state.dist_weights,theta_sins are up to date with x,y and thm
	//    respectively

	// Initialize random number generation
	unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
	std::default_random_engine generator(seed);
	angle_distribution<std::normal_distribution<double>> noise_dist(0, sqrt(2*state.T*state.dt));
	state.noise_sum = 0.0;

	double curnoise;
	for (int i = 0; i < state.conf.N; i++)
	{
		curnoise = noise_dist(generator);
		state.noise_sum += curnoise*curnoise;
		state.conf.thp[i] = angle_mod(state.conf.thm[i]
					    - state.dt*state.J * state.comp_Fterm(i)
				    	    + curnoise);
	}
}

void LeapFrogIntegrator::step_positions(State& state)
{
	//std::cout << "v0 = " << state.v0 << std::endl
		  //<< "dt = " << state.dt << std::endl
		  //<< "v0*dt = " << state.dt * state.v0 << std::endl
		  //<< "box = " << state.conf.box << std::endl
		  //<< "th0 = " << state.conf.thp[0] << std::endl
		  //<< std::endl;
	for(int i = 0; i < state.conf.N; i++)
	{
		state.conf.x[i] = mod_positive(state.conf.x[i] + state.v0*state.dt*cos(state.conf.thp[i]),
					       state.conf.box);
		state.conf.y[i] = mod_positive(state.conf.y[i] + state.v0*state.dt*sin(state.conf.thp[i]),
					       state.conf.box);
	}
}

void LeapFrogIntegrator::step(State& state)
{
	DBGOUT("stepping positions...");
	step_positions(state);
	state.update_pos_distances();

	double *tmp = state.conf.thm;
	state.conf.thm = state.conf.thp;
	state.conf.thp = tmp;
	tmp = 0;

	state.t += state.dt;

	DBGOUT("stepping orientations...");
	step_orientations(state);
	state.update_theta_distances();
}
