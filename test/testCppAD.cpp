// =============================================================================
//  CADET-semi-analytic - The semi analytic extension of
//  		CADET - The Chromatography Analysis and Design Toolkit
//  
//  Copyright © 2015-2017: Samuel Leweke¹
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
#include <vector>
#include <string>
#include <exception>
#include <algorithm>
#include <fstream>
#include <cmath>

#include "MPReal.hpp"
#include "MPComplex.hpp"
#include "MPAD.hpp"

#include <cppad/cppad.hpp>

#ifdef CASEMA_USE_FADBAD
	#include <tadiff.h>
#endif

template <typename real_t>
real_t func(const real_t& x)
{
	return real_t(1.0) / tanh(x);
}

template <typename real_t>
real_t func2(const real_t& x)
{
	return tanh(x);
}


#ifdef CASEMA_USE_FADBAD

	template <typename real_t>
	void checkFADBAD(const real_t& evalPos, std::size_t order, std::size_t digits)
	{
		fadbad::T<real_t> x;
		x[0] = evalPos;
		x[1] = 1.0;

		fadbad::T<mpfr::mpreal> val = func(x);
		val.eval(order);

		std::cout.flags(std::ios::scientific);
		std::cout.precision(digits);

		unsigned int fac = 1;

		std::cout << "========= FADBAD++ ==========" << std::endl;
		std::cout << "Taylor expansion of coth() at " << x[0] << " yields " << val[0] << ":\n";
		for (int i = 0; i < order; ++i)
		{
			if (i > 1)
				fac *= i;
			std::cout << "d^" << i << " f / dx^" << i << " = " << val[i] * fac << std::endl;;
		}
	}

#endif

template <typename real_t>
void checkCppAD(const real_t& evalPos, std::size_t order, std::size_t digits)
{
	std::vector< CppAD::AD<real_t> > x(1);
	x[0] = evalPos;

	// Start recording
	CppAD::Independent(x);
	// Execute function
	std::vector< CppAD::AD<real_t> > y(1);
	y[0] = func(x[0]);
	// Stop recording
	CppAD::ADFun<real_t> f(x, y);

	// Where to evaluate
	std::vector<real_t> pos(order+1);
	pos[0] = evalPos;
	pos[1] = real_t(1.0);
	for (std::size_t i = 2; i < order+1; ++i)
		pos[i] = real_t(0.0);

	const std::vector<real_t> dfdx = f.Forward(order, pos);

	std::cout.flags(std::ios::scientific);
	std::cout.precision(digits);

	unsigned int fac = 1;

	std::cout << "========= CppAD ==========" << std::endl;
	std::cout << "Taylor expansion of coth() at " << pos[0] << " yields " << dfdx[0] << ":\n";
	for (int i = 0; i < order; ++i)
	{
		if (i > 1)
			fac *= i;
		std::cout << "d^" << i << " f / dx^" << i << " = " << dfdx[i] * fac << std::endl;;
	}
}


int main(int argc, char** argv)
{
	const int precision = 100;
	mpfr::mpreal::set_default_prec(mpfr::digits2bits(precision));

	typedef mpfr::mpreal real_t;
	const std::size_t order = 10;

	const mpfr::mpreal pos = 1;

	checkCppAD(pos, order, 30);

#ifdef CASEMA_USE_FADBAD
	checkFADBAD(pos, order, 30);
#endif
	
	::mpfr_free_cache();

	return 0;
}

