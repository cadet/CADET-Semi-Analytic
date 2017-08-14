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

#include "MPAD.hpp"

#include <fadiff.h>
#include <badiff.h>
#include <tadiff.h>

#include "MPReal.hpp"
#include "MPComplex.hpp"

template <typename real_t>
real_t func(const real_t& x)
{
	return tanh(x);
}

struct CheckFunc
{
	static const char* Name(){ return "T tanh(const T&)"; }
	static const int m=1;
	static const int n=1;
	
	template <typename T>
	std::vector<T> operator()(const std::vector<T>& v)
	{
		std::vector<T> out(n);
		out[0]=func(v[0]);
		return out;
	}
	template <typename T>
	static std::vector<T> point()
	{
		std::vector< T > in(m);
		in[0] = T(1.0);
		return in;
	}
};


template <typename C>
struct Fdiff // Use stack-based forward differentiation
{
	static const int m=C::m;
	static const int n=C::n+C::m*C::n;
	
	template <typename T>
	std::vector<T> operator()(const std::vector<T>& v)
	{
		int i,j;
		std::vector< fadbad::F<T,C::m> > out(C::n);
		std::vector< fadbad::F<T,C::m> > in(C::m);
		for(i=0;i<C::m;++i)
		{
			in[i]=v[i];
			in[i].diff(i);
		}
		out=C()(in);
		std::vector<T> retval(n);
		for(j=0;j<C::n;++j)
		{
			retval[j]=out[j].x();
			for(i=0;i<C::m;++i)
			{
				retval[C::n+j*C::m+i]=out[j].d(i);
			}
		}
		return retval;
	}
};

template <typename C>
struct Bdiff
{
	static const int m=C::m;
	static const int n=C::n+C::m*C::n;
	
	template <typename T>
	std::vector<T> operator()(const std::vector<T>& v)
	{
		int i,j;
		std::vector< fadbad::B<T> > out(C::n);
		std::vector< fadbad::B<T> > in(C::m);
		for(i=0;i<C::m;++i)
		{
			in[i]=v[i];
		}
		out=C()(in);
		for(j=0;j<C::n;++j)
		{
			out[j].diff(j,C::n);
		}
		std::vector<T> retval(n);
		for(j=0;j<C::n;++j)
		{
			retval[j]=out[j].x();
			for(i=0;i<C::m;++i)
			{
				retval[C::n+j*C::m+i]=in[i].d(j);
			}
		}
		return retval;
	}
};


int main(int argc, char** argv)
{
	const int precision = 100;
	mpfr::mpreal::set_default_prec(mpfr::digits2bits(precision));

	// Taylor expansion
	fadbad::T<mpfr::mpreal> x;
	x = CheckFunc::point<mpfr::mpreal>()[0];
	x[1] = 1.0;

	fadbad::T<mpfr::mpreal> val = func(x);
	val.eval(12);

	std::cout.flags(std::ios::scientific);
	std::cout.precision(30);

	std::cout << "Taylor expansion of coth() at " << x[0] << " yields " << val[0] << ":\n";
	for (int i = 0; i < 12; ++i)
	{
		std::cout << "(1 / " << i << "!) * (d^" << i << " f / dx^" << i << ") = " << val[i] << std::endl;;
	}

	// Check derivatives by applying forward and backward mode multiple times
	{
		Fdiff<CheckFunc> Fd;
		std::vector<mpfr::mpreal> v = Fd(CheckFunc::point<mpfr::mpreal>());
		std::cout << "1st fwd: " << v[v.size()-1] << std::endl;
	}
	{
		Fdiff< Fdiff<CheckFunc> > Fd;
		std::vector<mpfr::mpreal> v = Fd(CheckFunc::point<mpfr::mpreal>());
		std::cout << "2nd fwd: " << v[v.size()-1] * 0.5 << std::endl;
	}
	{
		Fdiff< Fdiff< Fdiff<CheckFunc> > > Fd;
		std::vector<mpfr::mpreal> v = Fd(CheckFunc::point<mpfr::mpreal>());
		std::cout << "3rd fwd: " << v[v.size()-1] / (2.0*3.0) << std::endl;
	}
	{
		Fdiff< Fdiff< Fdiff< Fdiff<CheckFunc> > > > Fd;
		std::vector<mpfr::mpreal> v = Fd(CheckFunc::point<mpfr::mpreal>());
		std::cout << "4th fwd: " << v[v.size()-1] / (2.0*3.0*4.0) << std::endl;
	}

	{
		Bdiff<CheckFunc> Fd;
		std::vector<mpfr::mpreal> v = Fd(CheckFunc::point<mpfr::mpreal>());
		std::cout << "1st bwd: " << v[v.size()-1] << std::endl;
	}
	{
		Bdiff< Bdiff<CheckFunc> > Fd;
		std::vector<mpfr::mpreal> v = Fd(CheckFunc::point<mpfr::mpreal>());
		std::cout << "2nd bwd: " << v[v.size()-1] * 0.5 << std::endl;
	}
	{
		Bdiff< Bdiff< Bdiff<CheckFunc> > > Fd;
		std::vector<mpfr::mpreal> v = Fd(CheckFunc::point<mpfr::mpreal>());
		std::cout << "3rd bwd: " << v[v.size()-1] / (2.0 * 3.0) << std::endl;
	}
	{
		Bdiff< Bdiff< Bdiff< Bdiff<CheckFunc> > > > Fd;
		std::vector<mpfr::mpreal> v = Fd(CheckFunc::point<mpfr::mpreal>());
		std::cout << "4th bwd: " << v[v.size()-1] / (2.0*3.0*4.0) << std::endl;
	}

	::mpfr_free_cache();

	return 0;
}

