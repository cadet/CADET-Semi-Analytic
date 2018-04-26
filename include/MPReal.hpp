// =============================================================================
//  CADET-semi-analytic - The semi analytic extension of
//  		CADET - The Chromatography Analysis and Design Toolkit
//  
//  Copyright © 2015-2018: Samuel Leweke¹
//                                      
//    ¹ Forschungszentrum Juelich GmbH, IBG-1, Juelich, Germany.
//  
//  All rights reserved. This program and the accompanying materials
//  are made available under the terms of the GNU Public License v3.0 (or, at
//  your option, any later version) which accompanies this distribution, and
//  is available at http://www.gnu.org/licenses/gpl.html
// =============================================================================

#ifndef CASEMA_MPREAL_HPP_
#define CASEMA_MPREAL_HPP_

#include <limits>
#include <mpreal.h>

template <typename real_t>
inline real_t machineEpsilon(const real_t& v)
{
	return std::numeric_limits<real_t>::epsilon(v);
}

template <>
inline double machineEpsilon(const double& v)
{
	union helper_t
	{
		long long i64;
		double d64;
	} s;

	s.d64 = v;
	s.i64++;
	return (s.i64 < 0 ? v - s.d64 : s.d64 - v);
}

namespace std
{
	using mpfr::abs;
}

#endif
