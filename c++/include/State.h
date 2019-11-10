/* This file is distributed under MIT license as part of the project flocksims.
 * See LICENSE file for details.
 * (c) Simon Grosse-Holz, 2019
 */

# ifndef STATE_H
# define STATE_H

# include <string>
# include "H5Cpp.h"

namespace flocksims {

// class BaseIntegrator;
// class BaseReporter;

class Conformation
{
	public:
		Conformation(int N);
		Conformation(const Conformation&);
		~Conformation();

		void 		initialize_randomly();
		std::string	as_string();

		const int N;
		const double box = 1.0;
		double *x, *y, *thm, *thp;
};

class State
{
	friend class LeapFrogIntegrator;
	friend class HDF5Reporter;
	
	public:
		// This is the simulation state, so it needs to know some parameters...
		const double T, J, d_int, v0 = 1.0;
		double dt;

		State(const Conformation& conf, double T,
						double J,
						double dt,
						double d_int,
						double t = 0.0);
		State(const State&);
		~State();

		void		update_pos_distances();
		void		update_theta_distances();
		std::string	as_string() {return conf.as_string();}

		double	calc_polarization() const;
		double	calc_dSdt_discrete() const;
		double	calc_dSdt_expected() const;

		double	comp_Fterm(int i,
				   double *dist_weights=0,
				   double *theta_sins=0,
				   int N=-1) const;

	private:
		double t;
		Conformation conf;

		double *dist_weights, *theta_sins;
		double noise_sum;

		void	dpos_matrix(const double *x, const double *y, double *out,
				    int N=-1, double box=-1) const;
		void	dtheta_matrix(const double *theta, double *out, int N=-1) const;
};

}

# endif
