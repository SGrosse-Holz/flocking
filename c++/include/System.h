/* This file is distributed under MIT license as part of the project flocksims.
 * See LICENSE file for details.
 * (c) Simon Grosse-Holz, 2019
 */

# ifndef SYSTEM_H
# define SYSTEM_H

# include <string>

# include "State.h"
# include "Integrator.h"
# include "Reporter.h"

namespace flocksims {

class System
{
	public:
		System(State initial_state, BaseIntegrator *integrator);

		void	step(int n);
		State&	get_state() {return state;}
		void	report(BaseReporter *reporter, reportMode toReport);
		void	print();

	private:
		State		state;
		BaseIntegrator	*integrator;
};

}

# endif
