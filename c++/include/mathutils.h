/* This file is distributed under MIT license.
 * See LICENSE file for details.
 * (c) Simon Grosse-Holz, 2019
 */

// Note: by using templates, we prevent conflicts when this file is present in
// multiple translation units.

# ifndef MATHUTILS_SGH
# define MATHUTILS_SGH

# include <cmath>
# include <random>

// from  https://stackoverflow.com/questions/1903954/is-there-a-standard-sign-function-signum-sgn-in-c-c
template<typename T> int sgn(T val)
{
	return (T(0) < val) - (val < T(0 ));
}

template<typename T=double>
T mod_positive(T val, T by)
{
	return val - by*floor(val/by);
}

// Utilities for dealing with angles
// Angles are always in [-pi, pi).
template<typename T=double>
T angle_mod(T angle)
{
	return mod_positive(angle + M_PI, 2*M_PI) - M_PI;
}

template<typename T=double>
T angle_diff(T a1, T a2)
{
	return fmod(a1 - a2 + 3*M_PI, 2*M_PI) - M_PI;
}

template<typename T=double>
T angle_mean(T a1, T a2)
{
	T naive_diff = a1 - a2;
	T sum;
	if (naive_diff < -M_PI || naive_diff >= M_PI)
		sum = a1 + a2 + 2*M_PI;
	else
		sum = a1 + a2;

	return fmod(sum/2 + M_PI, 2*M_PI) - M_PI;
}

template<class base_dist>
class angle_distribution : public base_dist
{
	public:
		template<class... Args>
		angle_distribution(Args&&... args) : base_dist(args...) {}

		template<typename URNG>
		typename base_dist::result_type operator()(URNG& urng)
		{
			typename base_dist::result_type ret;
			while (true)
			{
				ret = base_dist::operator()(urng);
				if (ret >= -M_PI && ret < M_PI)
					break;
			}
			return ret;
		}
};

# endif
