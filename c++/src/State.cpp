/* This file is distributed under MIT license as part of the project flocksims.
 * See LICENSE file for details.
 * (c) Simon Grosse-Holz, 2019
 */

# include "State.h"
# include "mathutils.h"

# include <chrono>
# include <cmath>
# include <iomanip>
# include <iostream>
# include <random>
# include <sstream>

using namespace flocksims;

Conformation::Conformation(int N) : N(N)
{
	x = new double[N];
	y = new double[N];
	thm = new double[N];
	thp = new double[N];
}

Conformation::Conformation(const Conformation& other) : Conformation(other.N)
{
	for(int i = 0; i < N; i++)
	{
		x[i] = other.x[i];
		y[i] = other.y[i];
		thm[i] = other.thm[i];
		thp[i] = other.thp[i];
	}
}

Conformation::~Conformation()
{
	delete[] x;
	delete[] y;
	delete[] thm;
	delete[] thp;
}

void Conformation::initialize_randomly()
{
	unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
	std::default_random_engine generator(seed);
	std::uniform_real_distribution<double> pos_dist(0.0, box);
	std::uniform_real_distribution<double> ang_dist(-M_PI, M_PI);
	for(int i = 0; i < N; i++)
	{
		x[i] = pos_dist(generator);
		y[i] = pos_dist(generator);
		thm[i] = 0;
		thp[i] = ang_dist(generator);
	}
}

std::string Conformation::as_string()
{
	std::stringstream out = std::stringstream();
	out << "system size: " << N << std::endl;
	out << "theta-\tx\ty\ttheta+" << std::endl;
	for(int i = 0; i < N; i++)
	{
		out << thm[i] << "\t" << x[i] << "\t"
		    << y[i] << "\t" << thp[i] << std::endl;
	}
	return out.str();
}
	

State::State(const Conformation& conf, double T,
	       			       double J,
			       	       double dt,
				       double d_int,
				       double t)
	: T(T), J(J), dt(dt), d_int(d_int), t(t), conf(conf), noise_sum(0)
{
	dist_weights = new double[conf.N*conf.N];
	theta_sins = new double[conf.N*conf.N];

	update_pos_distances();
	update_theta_distances();
}

State::State(const State& other) : State(other.conf, other.T, other.J, other.dt, other.d_int, other.t)
{
	noise_sum = other.noise_sum;

	for (int i = 0; i < conf.N*conf.N; i++)
	{
		dist_weights[i] = other.dist_weights[i];
		theta_sins[i] = other.theta_sins[i];
	}
}

State::~State()
{
	delete[] dist_weights;
	delete[] theta_sins;
}

void State::dpos_matrix(const double *x, const double *y, double *out) const
{
	double xtmp, ytmp;
       	double halfbox = conf.box/2.0;
	int il, iu; // linear indices for lower and upper triangular matrix, respectively
	for(int i = 1; i < conf.N; i++)
	{
		for(int j = 0; j < i; j++)
		{
			il = i*conf.N + j;
			iu = j*conf.N + i;

			// Get PBC-aware distance and corresponding weights
			xtmp = fmod(x[i]-x[j] + halfbox, conf.box) - halfbox;
			ytmp = fmod(y[i]-y[j] + halfbox, conf.box) - halfbox;
			out[il] = sqrt(pow(xtmp, 2) + pow(ytmp, 2));
			out[iu] = exp(-out[il]/d_int);
		}
	}
}

void State::dtheta_matrix(const double *theta, double *out) const
{
	int il, iu;
	for (int i = 1; i < conf.N; i++)
	{
		for (int j = 0; j < i; j++)
		{
			il = i*conf.N + j;
			iu = j*conf.N + i;
			out[il] = theta[i] - theta[j]; // Attention: not in [-pi, pi)! Only use in periodic functions
			out[iu] = sin(out[il]);
		}
	}
}
	
void State::update_pos_distances()
{
	dpos_matrix(conf.x, conf.y, dist_weights);
}

void State::update_theta_distances()
{
	dtheta_matrix(conf.thp, theta_sins);
}

double State::calc_polarization() const
{
	double px=0, py=0;
	for (int i = 0; i < conf.N; i++)
	{
		px += cos(conf.thp[i]);
		py += sin(conf.thp[i]);
	}
	return sqrt((px*px + py*py))/conf.N;
}

double State::calc_dSdt_discrete() const
{
	double noise_sum_R = 0.0;
	double eta;
	int iu;
	for (int i = 0; i < conf.N; i++)
	{
		eta = 0.0;
		for (int j = 0; j < i; j++)
		{
			iu = j*conf.N + i;
			eta += dist_weights[iu]*theta_sins[iu];
		}
		for (int j = i+1; j < conf.N; j++)
		{
			iu = i*conf.N + j;
			eta -= dist_weights[iu]*theta_sins[iu]; // taking 'wrong' sin, therefore swap sign
		}
		eta = angle_mod(conf.thm - conf.thp + dt*J*eta);
		noise_sum_R += eta*eta;
	}

	return (noise_sum_R - noise_sum) / (4*T*dt);
}

double State::calc_dSdt_expected() const
{
	double *avg_theta_sins = new double[conf.N*conf.N];
	double *avg_theta = new double[conf.N];
	for (int i = 0; i < conf.N; i++)
		avg_theta[i] = angle_mean(conf.thm[i], conf.thp[i]);
	dtheta_matrix(avg_theta, avg_theta_sins);

	delete[] avg_theta;
	avg_theta = 0;

	double dS = 0.0;
	double factorj, tmpksum;
	int iuj, iuk;
	for (int i = 0; i < conf.N; i++)
	{
		// Always use k < j.
		// k < j < i
		// The sum over k is independent from j (except for the upper
		// bound), so we can compute and store it
		tmpksum = 0.0;
		factorj = avg_theta_sins[i]*dist_weights[i];
		dS += factorj*factorj;
		for (int j = 1; j < i; j++)
		{
			iuj = j*conf.N + i;
			factorj = avg_theta_sins[iuj]*dist_weights[iuj];
			iuk = (j-1)*conf.N + i;
			tmpksum += 2*avg_theta_sins[iuk]*dist_weights[iuk];
			dS += factorj*(factorj + tmpksum);
		}

		// i < j
		// The sum over k up to i is always the same, after that we do
		// the same as above
		tmpksum = 0.0;
		for (int k = 0; k < i; k++)
		{
			iuk = k*conf.N + i;
			tmpksum += 2*avg_theta_sins[iuk]*dist_weights[iuk];
		}
		iuj = i*conf.N + i+1;
		factorj = -avg_theta_sins[iuj]*dist_weights[iuj];
		dS += factorj*factorj;
		for (int j = i+2; j < conf.N; j++)
		{
			iuj = i*conf.N + j;
			factorj = -avg_theta_sins[iuj]*dist_weights[iuj]; // again 'wrong' sin
			iuk = i*conf.N + j-1;
			tmpksum -= 2*avg_theta_sins[iuk]*dist_weights[iuk]; // again 'wrong' sin
			dS += factorj*(factorj + tmpksum);
		}
	}

	delete[] avg_theta_sins;

	return J*J/T * dS/dt;
}
