/* This file is distributed under MIT license as part of the project flocksims.
 * See LICENSE file for details.
 * (c) Simon Grosse-Holz, 2019
 */

# ifndef INTEGRATOR_H
# define INTEGRATOR_H

# include "State.h"

namespace flocksims {

class Integrator
{
	public:
		Integrator(double T, double J, double d_int, double dt);

		void	step(State& state);
		void	step_positions(State& state);

	private:
		void	step_orientations(State& state);

		const double v0 = 1.0;
		const double noise_amp;
		const double J;
		const double d_int;
		const double dt;
};

}

# endif
