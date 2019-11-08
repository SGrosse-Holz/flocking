/* This file is distributed under MIT license as part of the project flocksims.
 * See LICENSE file for details.
 * (c) Simon Grosse-Holz, 2019
 */

# ifndef INTEGRATOR_H
# define INTEGRATOR_H

# include "State.h"

namespace flocksims {

class BaseIntegrator
{
	public:
		virtual ~BaseIntegrator() {}

		virtual void step(State& state) = 0;
};

class LeapFrogIntegrator : public BaseIntegrator
{
	public:
		void	step(State& state);

	private:
		void	step_orientations(State& state);
		void	step_positions(State& state);
};

}

# endif
