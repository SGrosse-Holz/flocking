# include "Integrator.h"

# include <chrono>
# include <cmath>
# include <functional>
# include <random>

using namespace flocksims;

Integrator::Integrator(double T, double J, double d_int, double dt) :
	noise_amp(sqrt(2*T*dt)), J(J), d_int(d_int), dt(dt)
{}

void Integrator::step_orientations(State& state)
{
	// TODO: optimize (maybe?)
	// Calculate all weights
	static double *dist_weights = new double[state.N*state.N]; // lower triangular: actual distances
								   // upper triangular: the corresponding weights
	static double *theta_diff = new double[state.N*state.N]; // similar
	static double *Dtheta_m = new double[state.N];
	static double entropy_sum = 0.0;

	state.dS_past = entropy_sum;
	double xtmp, ytmp, ntmp;
       	double halfbox = state.box/2.0;
	int lil, liu; // linear indices for lower and upper triangular matrix, respectively
	for(int i = 1; i < state.N; i++)
	{
		// For the backward probability, we combine the new D(theta)
		// and the old D(x,y).
		// Current status at this point:
		//  - state.dS_past = entropy_sum = (negative) sum of squares of the forward noise
		//  - Dtheta_m holds (minus) the last step in theta
		//  - dist_weights is still the same as for the previous step
		//  - theta_diff will be updated before the calculation
		ntmp = Dtheta_m[i];
		for(int j = 0; j < i; j++)
		{
			lil = i*state.N + j;
			liu = j*state.N + i;

			// Update angle differences
			theta_diff[lil] = state.theta[i]-state.theta[j];
			theta_diff[liu] = sin(theta_diff[lil]);

			// sum up the time-reversed noise
			ntmp += dist_weights[liu] * theta_diff[liu];

			// Get PBC-aware distance
			xtmp = abs(state.xp[i]-state.xp[j]);
			if(xtmp > halfbox) xtmp -= state.box;
			ytmp = abs(state.yp[i]-state.yp[j]);
			if(ytmp > halfbox) ytmp -= state.box;

			dist_weights[lil] = sqrt(pow(xtmp, 2)
					       + pow(ytmp, 2));

			dist_weights[liu] = J * exp(-dist_weights[lil]/d_int);
		}

		// At the end of the loop, this will be the entropy production
		// between this one and the last step (up to the prefactor).
		ntmp = fmod(ntmp+M_PI, 2*M_PI)-M_PI;
		state.dS_past += ntmp*ntmp;
	}
	state.dS_past /= 2*noise_amp*noise_amp;

	// Random initialization
	unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
	std::default_random_engine generator(seed);
	std::normal_distribution<double> noise_dist(0, noise_amp*noise_amp);
	auto noise = std::bind(noise_dist, generator);

	// Do the actual step
	// Since we have the distances separately, just overwrite in place
	entropy_sum = 0.0;
	double curnoise, theta_old;
	for(int i = 0; i < state.N; i++)
	{
		theta_old = state.theta[i];
		for(int j = 0; j < i; j++)
		{
			liu = j*state.N + i;
			state.theta[i] += -dt * dist_weights[liu] * theta_diff[liu];
		}
		for(int j = i+1; j < state.N; j++)
		{
			// Because of the triangular structure of the distance
			// matrix, here we have to access x_j - x_i instead of
			// x_i - x_j. For the spatial distance that's
			// irrelevant, because it is symmetric under this Z_2.
			// For theta the sign is important, so we get rid of
			// the global minus sign.
			liu = i*state.N + j;
			state.theta[i] += dt * dist_weights[liu] * theta_diff[liu];
		}
		curnoise = fmod(noise()+M_PI, 2*M_PI)-M_PI; // Get in [-M_PI, M_PI)
		entropy_sum -= curnoise*curnoise;
		state.theta[i] += curnoise;
		state.theta[i] = fmod(state.theta[i], 2*M_PI);
		Dtheta_m[i] = theta_old - state.theta[i]; // Doesn't need to be fmod'ed
	}
}

// Just fills in the state, i.e. calculates x+ from x- and theta
void Integrator::step_positions(State& state)
{
	for(int i = 0; i < state.N; i++)
	{
		state.xp[i] = fmod(state.xm[i] + v0*dt*cos(state.theta[i]), state.box);
		state.yp[i] = fmod(state.ym[i] + v0*dt*sin(state.theta[i]), state.box);
		if (state.xp[i] < 0) state.xp[i] += state.box;
		if (state.yp[i] < 0) state.yp[i] += state.box;
	}
}

// Assumes the state is consistent
// Accordingly, propagates orientations using x+
// 		does the timestep, i.e. x+ --> x-, t += dt
// 		calculates new x+
void Integrator::step(State& state)
{
	step_orientations(state);

	double *tmp = state.xm;
	state.xm = state.xp;
	state.xp = tmp;
	tmp = state.ym;
	state.ym = state.yp;
	state.yp = tmp;

	state.t += dt;

	step_positions(state);

	state.calculate_polarization();
}
