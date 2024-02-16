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

#include <vector>
#include <limits>
#include "MPReal.hpp"
#include "BesselZeros.hpp"

namespace
{
	inline bool equal(const mpfr::mpreal& ref, const mpfr::mpreal v)
	{
		return (abs(ref - v) <= std::numeric_limits<mpfr::mpreal>::epsilon()) || (abs(ref - v) / ref <= std::numeric_limits<mpfr::mpreal>::epsilon());
	}

	template <typename hankel_t>
	inline void dinisExpansionCoeffs(const std::vector<mpfr::mpreal>& zeros, std::vector<mpfr::mpreal>& coeff, hankel_t hankelCoeffs)
	{
		for (int i = 0; i < coeff.size(); ++i)
			coeff[i] = 2 * hankelCoeffs(i, zeros[i]) / sqr(mpfr::besselj0(zeros[i]));
	}

	inline mpfr::mpreal dinisExpansion(const std::vector<mpfr::mpreal>& zeros, const std::vector<mpfr::mpreal>& coeff, const mpfr::mpreal& p)
	{
		mpfr::mpreal r(0);
		for (int i = 0; i < zeros.size(); ++i)
			r += coeff[i] * mpfr::besselj0(p * zeros[i]);
		return r;
	}
}

TEST_CASE("Dini's expansion (r^2)", "[Bessel]")
{
	mpfr::mpreal::set_default_prec(mpfr::digits2bits(64));

	std::vector<mpfr::mpreal> zeros(10001);
	casema::besselZerosJ1(zeros.size(), zeros.data());

	std::vector<mpfr::mpreal> coeffs(zeros.size());
	dinisExpansionCoeffs(zeros, coeffs, [](int i, const mpfr::mpreal& z) -> mpfr::mpreal { return (i == 0) ? 0.25 : ( 2 * mpfr::besseljn(2, z) - z * mpfr::besseljn(3, z) ) / sqr(z); });

	std::vector<mpfr::mpreal> pts { 0.0, 0.125, 0.25, 0.375, 0.5, 0.625, 0.75, 0.875, 1.0 };

	for (const auto& p : pts)
	{
		CAPTURE(p);
		const mpfr::mpreal res = dinisExpansion(zeros, coeffs, p);
		const mpfr::mpreal ref = sqr(p);
		CAPTURE(res);
		CAPTURE(ref);
		CAPTURE(res-ref);

		CHECK(equal(res, ref));
	}
}

TEST_CASE("Dini's expansion (1-r^2)^-1", "[Bessel]")
{
	mpfr::mpreal::set_default_prec(mpfr::digits2bits(64));

	std::vector<mpfr::mpreal> zeros(10001);
	casema::besselZerosJ1(zeros.size(), zeros.data());

	std::vector<mpfr::mpreal> coeffs(zeros.size());
	dinisExpansionCoeffs(zeros, coeffs, [](int i, const mpfr::mpreal& z) -> mpfr::mpreal { return (i == 0) ? 1 : sin(z) / z; });

	std::vector<mpfr::mpreal> pts { 0.125, 0.25, 0.375, 0.5, 0.625, 0.75, 0.875 };

	for (const auto& p : pts)
	{
		CAPTURE(p);
		const mpfr::mpreal res = dinisExpansion(zeros, coeffs, p);
		const mpfr::mpreal ref = 1 / sqrt(1 - sqr(p));
		CAPTURE(res);
		CAPTURE(ref);
		CAPTURE(res-ref);

		CHECK(equal(res, ref));
	}
}
