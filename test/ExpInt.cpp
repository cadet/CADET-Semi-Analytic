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

#include <catch.hpp>

#include <limits>
#include "MPReal.hpp"
#include "ExpInt.hpp"

namespace
{
	inline bool equal(const mpfr::mpreal& ref, const mpfr::mpreal v)
	{
		return (abs(ref - v) <= std::numeric_limits<mpfr::mpreal>::epsilon() * 10) || (abs(ref - v) / ref <= std::numeric_limits<mpfr::mpreal>::epsilon() * 10);
	}
}

TEST_CASE("ExpInt and inverseExpInt", "[ExpInt]")
{
	// Test with 16, 32, 64 digits precision
	int prec = 16;
	for (int i = 0; i < 3; ++i)
	{
		mpfr::mpreal::set_default_prec(mpfr::digits2bits(prec));

		mpfr::mpreal x(2);
		for (int j = 0; j < 100; ++j)
		{
			const mpfr::mpreal y = casema::expInt(x);
			const mpfr::mpreal z = casema::inverseExpInt(y);

			CAPTURE(x);
			CAPTURE(prec);
			CAPTURE(z);
			CHECK(equal(z, x));

			x += 1;
		}

		x = "0.1";
		for (int j = 0; j < 5; ++j)
		{
			const mpfr::mpreal y = casema::expInt(x);
			const mpfr::mpreal z = casema::inverseExpInt(y);

			CAPTURE(x);
			CAPTURE(prec);
			CAPTURE(z);
			CHECK(equal(z, x));

			x /= 10;
		}

		x = 2;
		for (int j = 0; j < 100; ++j)
		{
			const mpfr::mpreal y = casema::inverseExpInt(x);
			const mpfr::mpreal z = casema::expInt(y);

			CAPTURE(x);
			CAPTURE(prec);
			CAPTURE(z);
			CHECK(equal(z, x));

			x += 1;
		}

		x = "0.1";
		for (int j = 0; j < 5; ++j)
		{
			const mpfr::mpreal y = casema::inverseExpInt(x);
			const mpfr::mpreal z = casema::expInt(y);

			CAPTURE(x);
			CAPTURE(prec);
			CAPTURE(z);
			CHECK(equal(z, x));

			x /= 10;
		}

		prec *= 2;
	}
}
