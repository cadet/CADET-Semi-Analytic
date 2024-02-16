// =============================================================================
//  CADET-semi-analytic - The semi-analytic extension of CADET
//  
//  Copyright © 2015-2020: Samuel Leweke¹²
//                                      
//    ¹ Forschungszentrum Juelich GmbH, IBG-1, Juelich, Germany.
//    ² University of Cologne, Cologne, Germany.
//  
//  All rights reserved. This program and the accompanying materials
//  are made available under the terms of the GNU Public License v3.0 (or, at
//  your option, any later version) which accompanies this distribution, and
//  is available at http://www.gnu.org/licenses/gpl.html
// =============================================================================

#ifndef CASEMA_CONSTANTS_HPP_
#define CASEMA_CONSTANTS_HPP_

#include <mpreal.h>

namespace casema
{

	template <typename T>
	struct Constants
	{
	};

	template <>
	struct Constants<double>
	{
		static double pi() { return 3.1415926535897932384626434; }
		static double eulerGamma() { return 0.5772156649015328606065121; }
		static void init(std::size_t precision) { }
		static void clear() { }
	};

	template <>
	struct Constants<mpfr::mpreal>
	{
		static mpfr::mpreal pi() { return mpfr::const_pi(); }
		static mpfr::mpreal eulerGamma() { return mpfr::const_euler(); }
		static void init(std::size_t precision) { mpfr::mpreal::set_default_prec(mpfr::digits2bits(precision)); }
		static void clear() { ::mpfr_free_cache(); }
	};

}

#endif
