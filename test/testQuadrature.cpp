// =============================================================================
//  CADET-semi-analytic - The semi analytic extension of
//  		CADET - The Chromatography Analysis and Design Toolkit
//  
//  Copyright © 2015-2019: Samuel Leweke¹
//                                      
//    ¹ Forschungszentrum Juelich GmbH, IBG-1, Juelich, Germany.
//  
//  All rights reserved. This program and the accompanying materials
//  are made available under the terms of the GNU Public License v3.0 (or, at
//  your option, any later version) which accompanies this distribution, and
//  is available at http://www.gnu.org/licenses/gpl.html
// =============================================================================

#include "CaSeMaConfig.hpp"

#include <iostream>

#include "MPReal.hpp"
#include <quadpack/workspace.hpp>


class PrecisionGuard
{
public:
	PrecisionGuard(std::ostream& os, std::size_t precision) : _os(os)
	{
		_curFlags = os.flags();
		_curPrec = os.precision();

		os.flags(std::ios::scientific);
		os.precision(precision);
	}

	~PrecisionGuard()
	{
		_os.precision(_curPrec);
		_os.flags(_curFlags);
	}

private:
	std::ostream& _os;
	std::ios_base::fmtflags _curFlags;
	std::streamsize _curPrec;
};


template <typename real_t>
struct Sin
{
	static const char* name() { return "Sin(x) on [0, 2*pi]"; }
	static real_t solution() { return real_t(0); }

	static real_t min() { return 0; }
	static real_t max() { return 2 * mpfr::const_pi(); }

	real_t operator()(const real_t& x) const
	{
		return sin(x);
	}
};


template <typename real_t>
struct HalfPi
{
	static const char* name() { return "1 / ((1-x)^2 + x^2)"; }
	static real_t solution() { return mpfr::const_pi() / real_t(2); }

	static real_t min() { return 0; }
	static real_t max() { return 1; }

	real_t operator()(const real_t& x) const
	{
		const real_t one(1);
		return one / (sqr(1-x) + sqr(x));
	}
};


template <typename real_t, template <class T> class Func>
void test(QuadPack::Workspace<real_t>& ws, const real_t& epsAbs, const real_t& epsRel)
{
	real_t result(0);
	real_t errorEst(0);
	const real_t sol = Func<real_t>::solution();

	Func<real_t> f;
	const int errorCode = ws.qag(f, Func<real_t>::min(), Func<real_t>::max(), epsAbs, epsRel, result, errorEst);

	std::cout << "========= " << Func<real_t>::name() << "=========" << std::endl;
	std::cout << "Return code: " << errorCode << " = " << QuadPack::returnCodeToString(errorCode) << std::endl;
	std::cout << "Result:   " << result << std::endl;
	std::cout << "Analytic: " << sol << std::endl;
	std::cout << "Est. error: " << errorEst << std::endl;
	std::cout << "Error:      " << abs(sol - result) << std::endl;
}


int main(int argc, char** argv)
{
	const int precision = 100;
	mpfr::mpreal::set_default_prec(mpfr::digits2bits(precision));

	typedef mpfr::mpreal real_t;

	const real_t epsAbs("1e-60");
	const real_t epsRel("1e-60");

	const std::size_t order = 10;
	const std::size_t maxIter = 50;

	std::cout << "Construct Gauss-Kronrod rule order " << order << std::endl;
	PrecisionGuard guard(std::cout, precision);

	QuadPack::Workspace<real_t> ws(maxIter, order);

	test<real_t, HalfPi>(ws, epsAbs, epsRel);
	test<real_t, Sin>(ws, epsAbs, epsRel);

	::mpfr_free_cache();

	return 0;
}

