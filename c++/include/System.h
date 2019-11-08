/* This file is distributed under MIT license as part of the project flocksims.
 * See LICENSE file for details.
 * (c) Simon Grosse-Holz, 2019
 */

# ifndef SYSTEM_H
# define SYSTEM_H

# include <string>

# include "State.h"
# include "Integrator.h"

namespace flocksims {

class System
{
	public:
				System(State initial_state, Integrator integrator);
		void		step(int n);
		std::string	state_string();
		const State&	get_state() {return state;}
		void		print();

	private:
		State		state;
		Integrator	integrator;
};

}

# endif
