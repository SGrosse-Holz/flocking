/* This file is distributed under MIT license as part of the project flocksims.
 * See LICENSE file for details.
 * (c) Simon Grosse-Holz, 2019
 */

# include "System.h"
# include <iostream>

using namespace flocksims;

System::System(State initial_state, BaseIntegrator *integrator)
	: state(initial_state), integrator(integrator)
{}

void System::step(int n)
{
	for (int i = 0; i < n; i++)
	{
		integrator->step(state);
	}
}

void System::print()
{
	std::cout << state.as_string() << std::endl;
}

void System::report(BaseReporter *reporter, reportMode toReport)
{
	reporter->report(state, toReport);
}
