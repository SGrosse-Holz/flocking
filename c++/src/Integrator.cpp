/* This file is distributed under MIT license as part of the project flocksims.
 * See LICENSE file for details.
 * (c) Simon Grosse-Holz, 2019
 */

# include "Integrator.h"
# include "mathutils.h"

# include <chrono>
# include <cmath>
# include <functional>
# include <random>

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
	auto noise = std::bind(noise_dist, generator);
	state.noise_sum = 0.0;

	int iu;
	double curnoise, *thp = state.conf.thp;
	for (int i = 0; i < state.conf.N; i++)
	{
		thp[i] = 0.0;
		for (int j = 0; j < i; j++)
		{
			iu = j*state.conf.N + i;
			thp[i] += state.dist_weights[iu]*state.theta_sins[iu];
		}
		for (int j = i+1; j < state.conf.N; j++)
		{
			iu = i*state.conf.N + j;
			thp[i] -= state.dist_weights[iu]*state.theta_sins[iu]; // c.f. State.cpp for sign
		}

		curnoise = noise();
		state.noise_sum += curnoise*curnoise;
		thp[i] = angle_mod(state.conf.thm[i] - state.dt*state.J*thp[i] + curnoise);
	}
}

void LeapFrogIntegrator::step_positions(State& state)
{
	for(int i = 0; i < state.conf.N; i++)
	{
		state.conf.x[i] = mod_positive(state.conf.x[i] + state.v0*state.dt*cos(state.conf.thp[i]), state.conf.box);
		state.conf.y[i] = mod_positive(state.conf.y[i] + state.v0*state.dt*sin(state.conf.thp[i]), state.conf.box);
	}
}

void LeapFrogIntegrator::step(State& state)
{
	step_positions(state);

	double *tmp = state.conf.thm;
	state.conf.thm = state.conf.thp;
	state.conf.thp = tmp;
	tmp = 0;

	state.t += state.dt;

	state.update_pos_distances();
	step_orientations(state);
	state.update_theta_distances();
}
