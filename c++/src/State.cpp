# include "State.h"

# include <chrono>
# include <cmath>
# include <iomanip>
# include <iostream>
# include <random>
# include <sstream>

using namespace flocksims;

State::State(int N, double t) : N(N), t(t)
{
	xp = new double[N];
	yp = new double[N];
	xm = new double[N];
	ym = new double[N];
	theta = new double[N];
}

State::State(const State& other) : State(other.N, other.t)
{
	polarization = other.polarization;
	dS_past = other.dS_past;

	for(int i = 0; i < N; i++)
	{
		xp[i] = other.xp[i];
		yp[i] = other.yp[i];
		xm[i] = other.xm[i];
		ym[i] = other.ym[i];
		theta[i] = other.theta[i];
	}
}

State::~State()
{
	delete[] xp;
	delete[] yp;
	delete[] xm;
	delete[] ym;
	delete[] theta;
}

// Note: initializes xp, yp, forward_noise to 0, those have to be set by the integrator!
void State::initialize_randomly()
{
	unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
	std::default_random_engine generator(seed);
	std::uniform_real_distribution<double> pos_dist(0.0, box);
	std::uniform_real_distribution<double> ang_dist(0.0, 2*M_PI);
	for(int i = 0; i < N; i++)
	{
		xp[i] = 0;
		yp[i] = 0;
		xm[i] = pos_dist(generator);
		ym[i] = pos_dist(generator);
		theta[i] = ang_dist(generator);
	}
}

std::string State::as_string()
{
	std::stringstream out = std::stringstream();
	out << "system size: " << N << ", t = " << std::setprecision(10) << t << std::endl;
	out << "x-\ty-\ttheta\tx+\ty+" << std::endl;
	for(int i = 0; i < N; i++)
	{
		out << xm[i] << "\t" << ym[i] << "\t"
		    << theta[i] << "\t"
		    << xp[i] << "\t" << yp[i] << std::endl;
	}
	return out.str();
}

void State::calculate_polarization()
{
	double px=0, py=0;
	for (int i = 0; i < N; i++)
	{
		px += cos(theta[i]);
		py += sin(theta[i]);
	}
	polarization = sqrt((px*px + py*py))/N;
}
