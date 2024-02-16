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

#include "BesselZeros.hpp"
#include "MPReal.hpp"

namespace
{
	inline mpfr::mpreal besselJ1derivative(const mpfr::mpreal& x)
	{
		return (mpfr::besselj0(x) - mpfr::besseljn(2, x)) / 2;
	}

	inline mpfr::mpreal besselZerosJ1newton(const mpfr::mpreal& guess, int maxIter = 1000)
	{
		mpfr::mpreal x = guess;
		for (int i = 0; i < maxIter; ++i)
		{
			const mpfr::mpreal res = mpfr::besselj1(x);
			const mpfr::mpreal dx = res / besselJ1derivative(x);

			if (abs(dx) <= machineEpsilon(x)) // (abs(res) <= std::numeric_limits<mpfr::mpreal>::epsilon())
				return x;

			x -= dx;
		}
		return x;
	}
}

namespace casema
{

void besselZerosJ1(int n, mpfr::mpreal* out) CASEMA_NOEXCEPT
{
	if (n <= 0)
		return;

	out[0] = 0;
	if (n == 1)
		return;

	const mpfr::mpreal pi = mpfr::const_pi();

	// Asymptotic expansion is J_1(x) = sqrt(2/(pi*x)) * cos(x - 3/4*pi)
	// Zeros of asymptotic expansion are near x = 5/4*pi + n * pi, n >= 0
	for (int i = 0; i < n-1; ++i)
	{
		// Use Newton iteration with initial guess from asymptotic expansion
		const mpfr::mpreal guess = (1.25 + i) * pi;
		out[i+1] = besselZerosJ1newton(guess);
	}
}

}
