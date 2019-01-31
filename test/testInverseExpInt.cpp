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
#include <limits>
#include "MPReal.hpp"

#include "ExpInt.hpp"

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
void writeLine(const real_t& val)
{
	const real_t in = casema::expInt(val);
	const real_t inv = casema::inverseExpInt(in);
	std::cout << val << "," << in << "," << inv << "," << abs(inv - val) << std::endl;
}


int main(int argc, char** argv)
{
	const int precision = 250;
	const int outPrec = 50;
	mpfr::mpreal::set_default_prec(mpfr::digits2bits(precision));

	typedef mpfr::mpreal real_t;
	PrecisionGuard guard(std::cout, outPrec);

	std::cout << "Eps = " << std::numeric_limits<real_t>::epsilon() << std::endl;
	std::cout << "Min = " << std::numeric_limits<real_t>::min() << std::endl;

	std::cout << "x        E_1       E_1^(-1)     error" << std::endl;

	if (argc <= 1)
	{
		const real_t x = 10;
		writeLine(x);
	}
	else
	{
		for (int i = 1; i < argc; ++i)
		{
			const real_t x(argv[i]);
			writeLine(x);
		}
	}

	::mpfr_free_cache();

	return 0;
}

