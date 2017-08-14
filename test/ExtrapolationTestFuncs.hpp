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

#ifndef CASEMA_TEST_EXTRAPOLATION_FUNCS_HPP_
#define CASEMA_TEST_EXTRAPOLATION_FUNCS_HPP_

#include "Constants.hpp"

template <typename real_t>
struct EulerE
{
	EulerE() : _n(0) { }

	const char* name() { return "Euler's number e"; }
	real_t solution() const { return exp(real_t(1)); }

	real_t next()
	{
		++_n;
		const real_t one(1);
		return pow(one  + one / _n, _n);
	}

	std::vector<real_t> sequence(std::size_t len) const
	{
		std::vector<real_t> seq(len, 0);
		const real_t one(1);
		for (std::size_t i = 0; i < seq.size(); ++i)
		{
			const real_t j = real_t(i);
			seq[i] = pow(one  + one / j, j);
		}
		return seq;
	}

	std::size_t _n;
};

template <typename real_t>
struct SquareRootOfTwo
{
	SquareRootOfTwo() : _cur(1) { }

	const char* name() { return "Square Root of Two"; }
	real_t solution() const { return sqrt(real_t(2)); }

	real_t next()
	{
		const real_t _prev = _cur;
		_cur = (_cur + real_t(2) / _cur) / real_t(2);
		return _prev;
	}

	std::vector<real_t> sequence(std::size_t len) const
	{
		std::vector<real_t> seq(len, 0);
		seq[0] = 1;
		for (std::size_t i = 0; i < seq.size()-1; ++i)
		{
			seq[i+1] = (seq[i] + real_t(2) / seq[i]) / real_t(2);
		}
		return seq;
	}

	real_t _cur;
};

template <typename real_t>
struct RootOfQuadraticPolynomial
{
	RootOfQuadraticPolynomial() : _cur(0) { }

	const char* name() { return "Root of x^2 - 4*x + 2"; }
	real_t solution() const { return real_t(2) - sqrt(real_t(2)); }

	real_t next()
	{
		const real_t _prev = _cur;
		_cur = (sqr(_cur) + real_t(2)) / real_t(4);
		return _prev;
	}

	std::vector<real_t> sequence(std::size_t len) const
	{
		std::vector<real_t> seq(len, 0);
		for (std::size_t i = 0; i < seq.size()-1; ++i)
		{
			seq[i+1] = (sqr(seq[i]) + real_t(2)) / real_t(4);
		}
		return seq;
	}

	real_t _cur;
};

template <typename real_t>
struct RootOfCubicPolynomial
{
	RootOfCubicPolynomial() : _cur("1.5") { }

	const char* name() { return "Root of x^3 + 4*x^2 - 10"; }
	real_t solution() const { return real_t(1) / real_t(3) * (real_t(-4) + cbrt(real_t(71) - real_t(3) * sqrt(real_t(105))) + cbrt(real_t(71) + real_t(3) * sqrt(real_t(105)))); }

	real_t next()
	{
		const real_t _prev = _cur;
		_cur = sqrt(real_t(10) - sqr(_cur) * _cur) / real_t(2);
		return _prev;
	}

	std::vector<real_t> sequence(std::size_t len) const
	{
		std::vector<real_t> seq(len, 0);
		seq[0] = real_t("1.5");
		for (std::size_t i = 0; i < seq.size()-1; ++i)
		{
			seq[i+1] = sqrt(real_t(10) - sqr(seq[i]) * seq[i]) / real_t(2);
		}
		return seq;
	}

	real_t _cur;
};

template <typename real_t>
struct RootOfExpEquation
{
	RootOfExpEquation() : _cur("0.3") { }

	const char* name() { return "Root of x - 3^x"; }
	real_t solution() const { return real_t("0.5478086216540974464505754081510218503459893377014890672937294550007263585900076805012650647619038191993101603827233709740722795994166359485635275881689621379811472888410464471969259623025787868324346387818769839657153664275805766648086828054855776490695149009539124841418137229456301630606785232340862402068698027572531351367335319227781321365012410951432459759405691961748386506830859539916401855541393799691716394262823235684721271696151326228809213229213977303216744465563068357955700718322080266602675860845555502729466197887490505045517473368191805662365178670551873342038788233915887316096947353571584722872096321991528721680224196810"); }

	real_t next()
	{
		const real_t _prev = _cur;
		_cur = exp(-_cur * log(real_t(3)));
		return _prev;
	}

	std::vector<real_t> sequence(std::size_t len) const
	{
		std::vector<real_t> seq(len, 0);
		seq[0] = real_t("0.3");
		for (std::size_t i = 0; i < seq.size()-1; ++i)
		{
			seq[i+1] = exp(-seq[i] * log(real_t(3)));
		}
		return seq;
	}

	real_t _cur;
};

template <typename real_t>
struct RiemannZeta
{
	RiemannZeta() : _cur(1), _n(0) { }

	const char* name() { return "Riemann Zeta at s = 2 (Pi^2 / 6)"; }
	real_t solution() const { return sqr(casema::Constants<real_t>::pi()) / real_t(6); }

	real_t next()
	{
		const real_t _prev = _cur;
		++_n;
		_cur = _cur + real_t(1) / sqr(real_t(_n+1));
		return _prev;
	}

	std::vector<real_t> sequence(std::size_t len) const
	{
		std::vector<real_t> seq(len, 0);
		seq[0] = real_t(1);
		for (std::size_t i = 0; i < seq.size()-1; ++i)
		{
			seq[i+1] = seq[i] + real_t(1) / sqr(real_t(i+2));
		}
		return seq;
	}

	real_t _cur;
	std::size_t _n;
};

template <typename real_t>
struct AlternatingHarmonicSeries
{
	AlternatingHarmonicSeries() : _cur(-1), _n(1) { }

	const char* name() { return "Alternating harmonic series = -ln(2)"; }
	real_t solution() const { return -log(real_t(2)); }

	real_t next()
	{
		const real_t _prev = _cur;
		++_n;
		if (_n % 2 == 0)
			_cur = _cur + real_t(1) / real_t(_n);
		else
			_cur = _cur - real_t(1) / real_t(_n);
		return _prev;
	}

	std::vector<real_t> sequence(std::size_t len) const
	{
		std::vector<real_t> seq(len, 0);
		seq[0] = -real_t(1);
		for (std::size_t i = 2; i <= seq.size(); ++i)
		{
			if (i % 2 == 0)
				seq[i-1] = seq[i-2] + real_t(1) / real_t(i);
			else
				seq[i-1] = seq[i-2] - real_t(1) / real_t(i);
		}
		return seq;
	}

	real_t _cur;
	std::size_t _n;
};

#endif
