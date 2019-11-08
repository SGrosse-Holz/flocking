/* This file is distributed under MIT license as part of the project flocksims.
 * See LICENSE file for details.
 * (c) Simon Grosse-Holz, 2019
 */

# ifndef STATE_H
# define STATE_H

# include <string>
# include "H5Cpp.h"

namespace flocksims {

class State {
	public:
		const int N;
		const double box = 1.0;

		double *xp, *yp, *xm, *ym, *theta; // These are gonna be int[N]
		double t;
		double polarization;
		double dS_past = 0.0;

		State(int N, double t = 0.0);
		State(const State&);
		~State();

		void 		initialize_randomly();
		std::string	as_string();
		void		append_to_hdf5_file(H5::H5File&, int ind=0);
		void		calculate_polarization();
};

}

# endif
