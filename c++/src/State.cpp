/* This file is distributed under MIT license as part of the project flocksims.
 * See LICENSE file for details.
 * (c) Simon Grosse-Holz, 2019
 */

# include "State.h"
# include "mathutils.h"

# include <chrono>
# include <cmath>
# include <fstream>
# include <iomanip>
# include <iostream>
# include <random>
# include <sstream>

using namespace flocksims;

Conformation::Conformation(int N, double density) : N(N), box(sqrt(N/density))
{
	x = new double[N];
	y = new double[N];
	thm = new double[N];
	thp = new double[N];
}

Conformation::Conformation(const Conformation& other) : Conformation(other.N, other.N/(other.box*other.box))
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

void Conformation::initialize_hbd(std::string file)
{
	std::ifstream ifile(file);
	if (!ifile.is_open())
		throw "Couldn't read file: "+file;

	std::string line;
	int i = 0;
	while (ifile >> line && i < N)
	{
		sscanf(line.c_str(), "%lf,%lf", x+i, y+i);
		i++;
	}

	if (i < N-1)
		throw "Didn't find enough positions in hbd file!";
	if (!ifile.eof())
		throw "There's positions left in the file!";

	unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
	std::default_random_engine generator(seed);
	std::uniform_real_distribution<double> ang_dist(-M_PI, M_PI);
	for (int j = 0; j < N; j++)
	{
		thm[j] = 0;
		thp[j] = ang_dist(generator);
	}
}

std::string Conformation::as_string()
{
	std::stringstream out = std::stringstream();
	out << "system size: " << N << " birds in a " << box << "x" << box << " box" << std::endl;
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
				       double v0,
				       double t)
	: T(T), J(J), d_int(d_int), v0(v0), dt(dt), t(t), conf(conf), noise_sum(0)
{
	dist_weights = new double[conf.N*conf.N];
	theta_sins = new double[conf.N*conf.N];

	update_pos_distances();
	update_theta_distances();
}

State::State(const State& other) : State(other.conf, other.T, other.J, other.dt, other.d_int, other.v0, other.t)
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

void State::dpos_matrix(const double *x, const double *y, double *out, int N, double box) const
{
	if (N < 0)
		N = conf.N;
	if (box < 0)
		box = conf.box;

	double xtmp, ytmp;
       	double halfbox = box/2.0;
	int il, iu; // linear indices for lower and upper triangular matrix, respectively
	for(int i = 1; i < N; i++)
	{
		for(int j = 0; j < i; j++)
		{
			il = i*N + j;
			iu = j*N + i;

			// Get PBC-aware distance and corresponding weights
			xtmp = mod_positive(x[i]-x[j] + halfbox, box) - halfbox;
			ytmp = mod_positive(y[i]-y[j] + halfbox, box) - halfbox;
			out[il] = sqrt(pow(xtmp, 2) + pow(ytmp, 2));
			out[iu] = exp(-out[il]/d_int);
		}
	}
}

void State::dtheta_matrix(const double *theta, double *out, int N) const
{
	if (N < 0)
		N = conf.N;

	int il, iu;
	for (int i = 1; i < N; i++)
	{
		for (int j = 0; j < i; j++)
		{
			il = i*N + j;
			iu = j*N + i;
			out[il] = theta[i] - theta[j]; // Attention: not in [-pi, pi)!
						       // Only use in periodic functions
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
	return sqrt(px*px + py*py)/conf.N;
}

double State::comp_Fterm(int i, double *dist_weights, double *theta_sins, int N) const
{
	if (!dist_weights)
		dist_weights = this->dist_weights;
	if (!theta_sins)
		theta_sins = this->theta_sins;
	if (N < 0)
		N = conf.N;

	double Fterm=0.0;
	int iu;

	Fterm = 0.0;
	for (int j = 0; j < i; j++)
	{
		iu = j*N + i;
		Fterm += dist_weights[iu]*theta_sins[iu];
	}
	for (int j = i+1; j < N; j++)
	{
		iu = i*N + j;
		Fterm -= dist_weights[iu]*theta_sins[iu]; // taking 'wrong' sin, therefore swap sign
	}

	return Fterm;
}

double State::calc_dSdt_discrete() const
{
	double eta;
	double noise_sum_R = 0.0;
	for (int i = 0; i < conf.N; i++)
	{
		eta = angle_mod(conf.thm[i] - conf.thp[i] + dt*J*comp_Fterm(i));
		noise_sum_R += eta*eta;
	}

	return (noise_sum_R - noise_sum) / (4*T*dt) / dt;
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

	double F=0.0, Fterm;
	for (int i = 0; i < conf.N; i++)
	{
		Fterm = comp_Fterm(i, dist_weights, avg_theta_sins);
		F += Fterm*Fterm;
	}

	delete[] avg_theta_sins;

	return J*J/T * F;
}
