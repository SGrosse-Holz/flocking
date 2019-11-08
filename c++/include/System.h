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
