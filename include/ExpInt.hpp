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

#ifndef CASEMA_EINT_HPP_
#define CASEMA_EINT_HPP_

#include <limits>
#include "Constants.hpp"

namespace casema
{

	template <typename real_t>
	real_t expInt(const real_t& x, unsigned int maxIter = 10000)
	{
		using std::abs;
//		using std::max;
		using std::exp;
		using std::log;

		if (x > 1)
		{
			// Implementation based on V. Pegoraro, P. Slusallek
			// "On the Evaluation of the Complex-Valued Exponential Integral", 
			// J. Graph. GPU, Game Tools. 15 (2011) 183–198. doi:10.1080/2151237X.2011.617177.

			// Use continued fraction expansion
			real_t c = real_t(0);
			real_t d = 1 + x;
			real_t e = exp(-x) / d;
			
			for (unsigned int k = 1; k <= maxIter; ++k)
			{
				const real_t common = 2 * k + 1 + x;
				c = real_t(1) / (common - k * k * c);
				d = common - k * k / d;
				const real_t delta = d * c;
				e /= delta;
				if (abs(delta - real_t(1)) <= std::numeric_limits<real_t>::epsilon())
				{
					return e;
				}
			}
		}
		else
		{
			// Implementation based on W.H. Press, S.A. Teukolsky, W.T. Vetterling, B.P. Flannery
			// "Numerical Recipes 3rd Edition: The Art of Scientific Computing", 3rd ed.,
			// Cambridge University Press, New York, NY, USA, 2007.

			// Use series expansion at x = 0
			real_t ans = -log(x) - Constants<real_t>::eulerGamma();
			real_t fact(1);

			for (unsigned int i = 1; i <= maxIter; ++i)
			{
				fact *= -x / i;
				const real_t del = -fact / i;
				ans += del;
	//            if (abs(del) <= std::numeric_limits<real_t>::epsilon())
				if (abs(del) <= abs(ans) * std::numeric_limits<real_t>::epsilon())
				{
					return ans;
				}
			}
		}

		// Failed to converge to desired relative accuracy
		return std::numeric_limits<real_t>::quiet_NaN();
	}


	template <typename real_t>
	real_t inverseExpInt(const real_t& x, const real_t& tol = std::numeric_limits<real_t>::epsilon(), const unsigned int maxIter = 10000)
	{
		using std::abs;
//		using std::max;
		using std::log;
		using std::expm1;
		using std::exp;

		// Apply Newton's method

		// Construct initial guess
		real_t p;
		if (x <= 1)
			p = max(real_t(1), -log(x));
		else
			p = real_t(1) / expm1(x);

		// Calculate residual
		real_t residual = casema::expInt(p, maxIter) - x;

		// Perform Newton steps until accuracy is reached
		for (unsigned int i = 0; (i < maxIter) && ((abs(residual) > tol) || (abs(residual) > tol * abs(x))); ++i)
		{
			const real_t step = exp(p) * residual;

			// Stop if adding delta does not yield any change in p
			if (abs(p * step) <= machineEpsilon(p))
				break;

			// Damping
			real_t alpha(1);

			if (alpha * step <= -1)
				alpha = -real_t(1) / step + tol;

			// Try line search with at most 10 steps to decrease residual
			real_t cand = p * (real_t(1) + alpha * step);
			real_t newResidual = casema::expInt(cand, maxIter) - x;
			for (unsigned int i = 0; (i < 10) && (abs(newResidual) > abs(residual)); ++i)
			{
				alpha *= real_t(0.5);
				cand = p * (real_t(1) + alpha * step);
				newResidual = casema::expInt(cand, maxIter) - x;
			}

			// Accept step
			p = cand;
			residual = newResidual;
		}

		return p;
	}

}

#endif
