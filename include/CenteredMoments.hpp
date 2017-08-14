// =============================================================================
//  CADET-semi-analytic - The semi analytic extension of
//  		CADET - The Chromatography Analysis and Design Toolkit
//  
//  Copyright © 2015: Samuel Leweke¹
//                                      
//    ¹ Forschungszentrum Juelich GmbH, IBG-1, Juelich, Germany.
//  
//  All rights reserved. This program and the accompanying materials
//  are made available under the terms of the GNU Public License v3.0 (or, at
//  your option, any later version) which accompanies this distribution, and
//  is available at http://www.gnu.org/licenses/gpl.html
// =============================================================================

#ifndef CASEMA_MPREAL_CENTMOM_HPP_
#define CASEMA_MPREAL_CENTMOM_HPP_

#include "MPReal.hpp"
#include <vector>

namespace mpfr
{
	inline mpreal binom(unsigned long int n, unsigned long int k)
	{
		mpz_t ret;
		mpz_init(ret);

		mpz_bin_uiui(ret, n, k);

		mpreal y(ret);
		mpz_clear(ret);

		return y;
	}
}

namespace casema
{

	template <typename T>
	std::vector<T> centralMomentsFromNonCentral(const std::vector<T>& nonCent)
	{
		std::vector<T> cent(nonCent.size(), T(0));

		// Compute 0th and 1st moment by copying
		for (std::size_t n = 0; n < std::min(static_cast<std::size_t>(2), nonCent.size()); ++n)
			cent[n] = nonCent[n];

		// Use binomial expansion for higher order moments
		for (std::size_t n = 2; n < nonCent.size(); ++n)
		{
			// Use Horner scheme for evaluating the polynomial in \mu_1

			// Compute highest order term (-1)^(n-1) * \mu_1^n * (n-1)
			if (n % 2 == 0)
				cent[n] = -nonCent[1] * T(n-1);
			else
				cent[n] = nonCent[1] * T(n-1);

			// Compute rest of the polynomial
			for (std::size_t k = n-2; k >= 0; --k)
			{
				if (k % 2 == 0)
					cent[n] = nonCent[1] * cent[n] + mpfr::binom(n,k) * nonCent[n-k];
				else
					cent[n] = nonCent[1] * cent[n] - mpfr::binom(n,k) * nonCent[n-k];

				if (k == 0)
					break;
			}
		}
		return cent;
	}

}

#endif
