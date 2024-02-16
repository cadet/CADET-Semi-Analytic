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

/**
 * @file 
 * Provides tools for managing scopes of IParameterProvider
 */

#ifndef CASEMA_PARAMREADERUTILS_HPP_
#define CASEMA_PARAMREADERUTILS_HPP_

#include "io/ParameterProvider.hpp"
#include "Exceptions.hpp"
#include "MPReal.hpp"

#include <vector>

namespace casema
{

	inline std::vector<mpfr::mpreal> toMPreal(const std::vector<double>& v)
	{
		std::vector<mpfr::mpreal> d(0);
		d.reserve(v.size());
		for (double vv : v)
			d.push_back(vv);

		return d;
	}

	inline std::vector<mpfr::mpreal> toMPreal(const std::vector<double>& v, int n)
	{
		std::vector<mpfr::mpreal> d(0);
		d.reserve(n);
		if (v.size() == 1)
		{
			for (int i = 0; i < n; ++i)
				d.push_back(v[0]);
		}
		else
		{
			for (double vv : v)
				d.push_back(vv);
		}

		return d;
	}

} // namespace casema

#endif  // CASEMA_PARAMREADERUTILS_HPP_
