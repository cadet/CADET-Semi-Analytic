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

#include "DurbinsMethod.hpp"

#include <vector>

TEST_CASE("Inverse Laplace exp(-2 * sqrt(s)) single", "[Durbin]")
{
	const std::size_t precision = 50;
	mpfr::mpreal::set_default_prec(mpfr::digits2bits(precision));

	const double tMaxDbl = 10.0;
	const mpfr::mpreal tMax = tMaxDbl;
	const mpfr::mpreal error("1e-25");
	const mpfr::mpreal weight("0.005");
/*
	const std::size_t precision = 100;
	mpfr::mpreal::set_default_prec(mpfr::digits2bits(precision));

	const double tMaxDbl = 10.0;
	const mpfr::mpreal tMax = tMaxDbl;
	const mpfr::mpreal error("1e-50");
	const mpfr::mpreal weight("0.005");
*/

	const mpfr::mpreal T = tMax / 2 * mpfr::mpreal(1.01);
	const mpfr::mpreal abscissaSafety = 0.25;
	const mpfr::mpreal abscissa = log1p(0.25 / (weight * error)) / (2 * T) + abscissaSafety;

	// exp(-sqrt(N)) * (1 + sqrt(N)) <= (1-weight) / 2 * error * T / exp(2 * abscissa * T)
	const mpfr::mpreal rhs = (1-weight) / 2 * error * T / exp(2 * abscissa * T);
	mpfr::mpreal nmp = sqr(log(rhs));
	for (int i = 0; i < 1000; ++i)
	{
		const mpfr::mpreal dx = 2 * (rhs * exp(sqrt(nmp)) - 1 - sqrt(nmp));
		if (abs(dx) <= machineEpsilon(nmp))
			break;

		nmp -= dx;
	}

	const std::size_t nSummands = static_cast<std::size_t>(ceil(nmp));

	std::vector<mpfr::mpreal> timeGrid(99);
	for (int i = 0; i < timeGrid.size(); ++i)
		timeGrid[i] = (i + 1.0) * tMaxDbl / timeGrid.size();

	const casema::util::SlicedVector<mpfr::mpreal> result = casema::invertLaplace(
		[](const mpfr::mpcomplex& s, mpfr::mpcomplex* r) { r[0] = exp(-2 * sqrt(s)); },
		1, nSummands, precision, tMax, abscissa, timeGrid.data(), timeGrid.size());

	for (int j = 0; j < timeGrid.size(); ++j)
	{
		const mpfr::mpreal& t = timeGrid[j];
		const mpfr::mpreal ref = exp(-1.0 / t) / sqrt(mpfr::const_pi() * sqr(t)*t);

//		std::cout << "t: " << t << " lap: " << result(0, j) << " ref: " << ref << " diff " << abs(ref - result(0, j)) << " rel " << abs(ref - result(0, j)) / abs(ref) << "\n";

		CAPTURE(j);
		CHECK(abs(ref - result(0, j)) <= 1e-24);
	}
}

TEST_CASE("Inverse Laplace exp(-2 * sqrt(s)) two outlets", "[Durbin]")
{
	const std::size_t precision = 50;
	mpfr::mpreal::set_default_prec(mpfr::digits2bits(precision));

	const double tMaxDbl = 10.0;
	const mpfr::mpreal tMax = tMaxDbl;
	const mpfr::mpreal error("1e-25");
	const mpfr::mpreal weight("0.005");
/*
	const std::size_t precision = 100;
	mpfr::mpreal::set_default_prec(mpfr::digits2bits(precision));

	const double tMaxDbl = 10.0;
	const mpfr::mpreal tMax = tMaxDbl;
	const mpfr::mpreal error("1e-50");
	const mpfr::mpreal weight("0.005");
*/

	const mpfr::mpreal T = tMax / 2 * mpfr::mpreal(1.01);
	const mpfr::mpreal abscissaSafety = 0.25;
	const mpfr::mpreal abscissa = log1p(0.25 / (weight * error)) / (2 * T) + abscissaSafety;

	// exp(-sqrt(N)) * (1 + sqrt(N)) <= (1-weight) / 2 * error * T / exp(2 * abscissa * T)
	const mpfr::mpreal rhs = (1-weight) / 2 * error * T / exp(2 * abscissa * T);
	mpfr::mpreal nmp = sqr(log(rhs));
	for (int i = 0; i < 1000; ++i)
	{
		const mpfr::mpreal dx = 2 * (rhs * exp(sqrt(nmp)) - 1 - sqrt(nmp));
		if (abs(dx) <= machineEpsilon(nmp))
			break;

		nmp -= dx;
	}

	const std::size_t nSummands = static_cast<std::size_t>(ceil(nmp));

	std::vector<mpfr::mpreal> timeGrid(99);
	for (int i = 0; i < timeGrid.size(); ++i)
		timeGrid[i] = (i + 1.0) * tMaxDbl / timeGrid.size();

	const casema::util::SlicedVector<mpfr::mpreal> result = casema::invertLaplace(
		[](const mpfr::mpcomplex& s, mpfr::mpcomplex* r)
		{
			const mpfr::mpcomplex v = exp(-2 * sqrt(s));
			r[0] = v;
			r[1] = -v;
		},
		2, nSummands, precision, tMax, abscissa, timeGrid.data(), timeGrid.size());

	for (int j = 0; j < timeGrid.size(); ++j)
	{
		const mpfr::mpreal& t = timeGrid[j];
		const mpfr::mpreal ref = exp(-1.0 / t) / sqrt(mpfr::const_pi() * sqr(t)*t);

//		std::cout << "t: " << t << " lap: " << result(0, j) << " ref: " << ref << " diff " << abs(ref - result(0, j)) << " rel " << abs(ref - result(0, j)) / abs(ref) << "\n";

		CAPTURE(j);
		CHECK(abs(ref - result(0, j)) <= 1e-24);
		CHECK(abs(ref + result(1, j)) <= 1e-24);
	}
}
