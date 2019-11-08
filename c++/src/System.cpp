# include "System.h"
# include <iostream>

using namespace flocksims;

System::System(State initial_state, Integrator integrator) : state(initial_state), integrator(integrator)
{
	integrator.step_positions(state);
}

void System::step(int n)
{
	for (int i = 0; i < n; i++)
	{
		integrator.step(state);
	}
}

void System::print()
{
	std::cout << state.as_string() << std::endl;
}
