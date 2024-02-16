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

#ifndef CASEMA_MPCOMPLEX_HPP_
#define CASEMA_MPCOMPLEX_HPP_

#include <string>
#include <ostream>

#include <mpc.h>
#include <mpreal.h>

namespace mpfr
{
	
	class mpcomplex
	{
	public:
		mpcomplex() 
		{ 
			::mpc_init2(_val, mpreal::get_default_prec());
			::mpc_set_d(_val, 0.0, MPC_RNDNN);
		}

		mpcomplex(const mpcomplex& v) 
		{ 
			::mpc_init2(_val, mpreal::get_default_prec());
			::mpc_set(_val, v._val, MPC_RNDNN);
		}

		mpcomplex(mpcomplex&& v) noexcept
		{
			mpfr_set_uninitialized(_val->re);
			mpc_swap(_val, v._val);
		}

		mpcomplex(const mpreal& v) 
		{ 
			::mpc_init2(_val, mpreal::get_default_prec());
			::mpc_set_fr(_val, v.mpfr_ptr(), MPC_RNDNN);
		}

		mpcomplex(int v)
		{
			::mpc_init2(_val, mpreal::get_default_prec());
			::mpc_set_si(_val, v, MPC_RNDNN);
		}

		mpcomplex(long int v)
		{
			::mpc_init2(_val, mpreal::get_default_prec());
			::mpc_set_si(_val, v, MPC_RNDNN);
		}

		mpcomplex(unsigned int v)
		{
			::mpc_init2(_val, mpreal::get_default_prec());
			::mpc_set_ui(_val, v, MPC_RNDNN);
		}

		mpcomplex(unsigned long int v)
		{
			::mpc_init2(_val, mpreal::get_default_prec());
			::mpc_set_ui(_val, v, MPC_RNDNN);
		}

		mpcomplex(double v)
		{
			::mpc_init2(_val, mpreal::get_default_prec());
			::mpc_set_d(_val, v, MPC_RNDNN);
		}

		mpcomplex(const mpreal& vr, const mpreal& vi) 
		{ 
			::mpc_init2(_val, mpreal::get_default_prec());
			::mpc_set_fr_fr(_val, vr.mpfr_ptr(), vi.mpfr_ptr(), MPC_RNDNN);
		}

		mpcomplex(double vr, double vi)
		{
			::mpc_init2(_val, mpreal::get_default_prec());
			::mpc_set_d_d(_val, vr, vi, MPC_RNDNN);
		}

		mpcomplex(const char* const str)
		{
			mpreal x(str);
			::mpc_init2(_val, mpreal::get_default_prec());
			::mpc_set_fr(_val, x.mpfr_ptr(), MPC_RNDNN);
		}

		mpcomplex(const std::string& str)
		{
			mpreal x(str);
			::mpc_init2(_val, mpreal::get_default_prec());
			::mpc_set_fr(_val, x.mpfr_ptr(), MPC_RNDNN);
		}

		~mpcomplex()
		{
			if (mpfr_is_initialized(_val->re))
				mpc_clear(_val);
		}

		mpreal real() const
		{ 
			mpreal x; 
			mpc_real(x.mpfr_ptr(), _val, mpreal::get_default_rnd());
			return x;
		}

		mpreal imag() const
		{ 
			mpreal x; 
			mpc_imag(x.mpfr_ptr(), _val, mpreal::get_default_rnd());
			return x;
		}

		friend const mpcomplex operator-(const mpcomplex& a);

		mpcomplex& operator=(const mpreal& v);
		mpcomplex& operator=(const mpcomplex& v);
		mpcomplex& operator=(const double v);        
		mpcomplex& operator=(const char* s);
		mpcomplex& operator=(const std::string& s);
		mpcomplex& operator=(mpcomplex&& v) noexcept;

		// +
		mpcomplex& operator+=(const mpcomplex& v);
		mpcomplex& operator+=(const mpreal& v);
		mpcomplex& operator+=(const double u);
		friend const mpcomplex operator+(const mpreal& b,     const mpcomplex& a);
		friend const mpcomplex operator+(const mpcomplex& b,  const mpreal& a);
		friend const mpcomplex operator+(const mpcomplex& b,  const mpcomplex& a);

		// -
		mpcomplex& operator-=(const mpcomplex& v);
		mpcomplex& operator-=(const mpreal& v);
		mpcomplex& operator-=(const double u);
		friend const mpcomplex operator-(const mpreal& b,     const mpcomplex& a);
		friend const mpcomplex operator-(const mpcomplex& b,  const mpreal& a);
		friend const mpcomplex operator-(const mpcomplex& b,  const mpcomplex& a);

		// *
		mpcomplex& operator*=(const mpcomplex& v);
		mpcomplex& operator*=(const mpreal& v);
		mpcomplex& operator*=(const double v);
		friend const mpcomplex operator*(const mpreal& b,     const mpcomplex& a);
		friend const mpcomplex operator*(const mpcomplex& b,  const mpreal& a);
		friend const mpcomplex operator*(const mpcomplex& b,  const mpcomplex& a);
		friend const mpcomplex operator*(const int b,         const mpcomplex& a);
		friend const mpcomplex operator*(const mpcomplex& b,  const int a);
		
		// /
		mpcomplex& operator/=(const mpcomplex& v);
		mpcomplex& operator/=(const mpreal& v);
		mpcomplex& operator/=(const double v);
		friend const mpcomplex operator/(const mpreal& b,      const mpcomplex& a);
		friend const mpcomplex operator/(const mpcomplex& b,   const mpreal& a);
		friend const mpcomplex operator/(const mpcomplex& b,   const mpcomplex& a);
		friend const mpcomplex operator/(const unsigned int b, const mpcomplex& a);
		friend const mpcomplex operator/(const mpcomplex& b,   const unsigned int a);

		// Boolean Operators
		friend bool operator== (const mpcomplex& a, const mpcomplex& b);
		friend bool operator!= (const mpcomplex& a, const mpcomplex& b);

		// Optimized specializations for boolean operators
		friend bool operator== (const mpcomplex& a, const long int b);

		// Convert mpcomplex to string with n significant digits in base b
		// n = -1 -> convert with the maximum available digits
		std::string toString(int n = -1, int b = 10) const;

		std::ostream& output(std::ostream& os) const;

		// Get raw pointers so that mpreal can be directly used in raw mpfr_* functions
		::mpc_ptr    mpc_ptr() { return _val; }
		::mpc_srcptr mpc_ptr()    const  { return _val; }
		::mpc_srcptr mpc_srcptr() const  { return _val; }

	protected:
		::mpc_t _val;
	};


	inline const mpcomplex operator-(const mpcomplex& a)
	{
		mpcomplex c(a);
		mpc_neg(c._val, a.mpc_srcptr(), MPC_RNDNN);
		return c;
	}

	//////////////////////////////////////////////////////////////////////////
	// = Assignment
	inline mpcomplex& mpcomplex::operator=(const mpcomplex& v)
	{
		if (this != &v)
		{
			mpc_set(_val, v._val, MPC_RNDNN);
		}
		return *this;
	}

	inline mpcomplex& mpcomplex::operator=(const mpreal& v)
	{
		mpc_set_fr(_val, v.mpfr_srcptr(), MPC_RNDNN);
		return *this;
	}

	inline mpcomplex& mpcomplex::operator=(const double v)
	{
		mpc_set_d(_val, v, MPC_RNDNN);
		return *this;
	}

	inline mpcomplex& mpcomplex::operator=(mpcomplex&& v) noexcept
	{
		mpc_swap(_val, v._val);
		return *this;
	}

	//////////////////////////////////////////////////////////////////////////
	// + Addition
	inline mpcomplex& mpcomplex::operator+=(const mpcomplex& v)
	{
		mpc_add(_val, _val, v._val, MPC_RNDNN);
		return *this;
	}

	inline mpcomplex& mpcomplex::operator+=(const mpreal& v)
	{
		mpc_add_fr(_val, _val, v.mpfr_srcptr(), MPC_RNDNN);
		return *this;
	}

	inline mpcomplex& mpcomplex::operator+= (const double u)
	{
		mpreal x(u);
		mpc_add_fr(_val,_val, x.mpfr_srcptr(),MPC_RNDNN);
		return *this;
	}


	inline const mpcomplex operator+(const mpcomplex& a, const mpreal& b)
	{
		mpcomplex c;
		mpc_add_fr(c._val, a._val, b.mpfr_srcptr(), MPC_RNDNN);
		return c;
	}

	inline const mpcomplex operator+(const mpreal& b, const mpcomplex& a)
	{
		mpcomplex c;
		mpc_add_fr(c._val, a._val, b.mpfr_srcptr(), MPC_RNDNN);
		return c;
	}

	inline const mpcomplex operator+(const mpcomplex& a, const mpcomplex& b)
	{
		mpcomplex c;
		mpc_add(c._val, a._val, b._val, MPC_RNDNN);
		return c;
	}

	//////////////////////////////////////////////////////////////////////////
	// - Subtraction
	inline mpcomplex& mpcomplex::operator-=(const mpcomplex& v)
	{
		mpc_sub(_val, _val, v._val, MPC_RNDNN);
		return *this;
	}

	inline mpcomplex& mpcomplex::operator-=(const mpreal& v)
	{
		mpc_sub_fr(_val, _val, v.mpfr_srcptr(), MPC_RNDNN);
		return *this;
	}

	inline mpcomplex& mpcomplex::operator-= (const double u)
	{
		mpreal x(u);
		mpc_sub_fr(_val,_val, x.mpfr_srcptr(),MPC_RNDNN);
		return *this;
	}


	inline const mpcomplex operator-(const mpcomplex& a, const mpreal& b)
	{
		mpcomplex c;
		mpc_sub_fr(c._val, a._val, b.mpfr_srcptr(), MPC_RNDNN);
		return c;
	}

	inline const mpcomplex operator-(const mpreal& a, const mpcomplex& b)
	{
		mpcomplex c;
		mpc_sub_fr(c._val, b._val, a.mpfr_srcptr(), MPC_RNDNN);
		mpc_neg(c._val, c._val, MPC_RNDNN);
		return c;
	}

	inline const mpcomplex operator-(const mpcomplex& a, const mpcomplex& b)
	{
		mpcomplex c;
		mpc_sub(c._val, a._val, b._val, MPC_RNDNN);
		return c;
	}

	//////////////////////////////////////////////////////////////////////////
	// * Multiplication
	inline mpcomplex& mpcomplex::operator*=(const mpcomplex& v)
	{
		mpc_mul(_val, _val, v._val, MPC_RNDNN);
		return *this;
	}

	inline mpcomplex& mpcomplex::operator*=(const mpreal& v)
	{
		mpc_mul_fr(_val, _val, v.mpfr_srcptr(), MPC_RNDNN);
		return *this;
	}

	inline mpcomplex& mpcomplex::operator*= (const double u)
	{
		mpreal x(u);
		mpc_mul_fr(_val,_val, x.mpfr_srcptr(),MPC_RNDNN);
		return *this;
	}


	inline const mpcomplex operator*(const mpcomplex& a, const mpreal& b)
	{
		mpcomplex c;
		mpc_mul_fr(c._val, a._val, b.mpfr_srcptr(), MPC_RNDNN);
		return c;
	}

	inline const mpcomplex operator*(const mpreal& a, const mpcomplex& b)
	{
		mpcomplex c;
		mpc_mul_fr(c._val, b._val, a.mpfr_srcptr(), MPC_RNDNN);
		return c;
	}

	inline const mpcomplex operator*(const mpcomplex& a, const mpcomplex& b)
	{
		mpcomplex c;
		mpc_mul(c._val, a._val, b._val, MPC_RNDNN);
		return c;
	}

	inline const mpcomplex operator*(const mpcomplex& a, const int b)
	{
		mpcomplex c;
		mpc_mul_si(c._val, a._val, b, MPC_RNDNN);
		return c;
	}

	inline const mpcomplex operator*(const int a, const mpcomplex& b)
	{
		mpcomplex c;
		mpc_mul_si(c._val, b._val, a, MPC_RNDNN);
		return c;
	}

	//////////////////////////////////////////////////////////////////////////
	// / Division
	inline mpcomplex& mpcomplex::operator/=(const mpcomplex& v)
	{
		mpc_div(_val, _val, v._val, MPC_RNDNN);
		return *this;
	}

	inline mpcomplex& mpcomplex::operator/=(const mpreal& v)
	{
		mpc_div_fr(_val, _val, v.mpfr_srcptr(), MPC_RNDNN);
		return *this;
	}

	inline mpcomplex& mpcomplex::operator/= (const double u)
	{
		mpreal x(u);
		mpc_div_fr(_val,_val, x.mpfr_srcptr(),MPC_RNDNN);
		return *this;
	}


	inline const mpcomplex operator/(const mpcomplex& a, const mpreal& b)
	{
		mpcomplex c;
		mpc_div_fr(c._val, a._val, b.mpfr_srcptr(), MPC_RNDNN);
		return c;
	}

	inline const mpcomplex operator/(const mpreal& a, const mpcomplex& b)
	{
		mpcomplex c(a);
		mpc_div(c._val, c._val, b._val, MPC_RNDNN);
		return c;
	}

	inline const mpcomplex operator/(const mpcomplex& a, const mpcomplex& b)
	{
		mpcomplex c;
		mpc_div(c._val, a._val, b._val, MPC_RNDNN);
		return c;
	}

	inline const mpcomplex operator/(const mpcomplex& a, const unsigned int b)
	{
		mpcomplex c;
		mpc_div_ui(c._val, a._val, b, MPC_RNDNN);
		return c;
	}

	inline const mpcomplex operator/(const unsigned int a, const mpcomplex& b)
	{
		mpcomplex c;
		mpc_ui_div(c._val, a, b._val, MPC_RNDNN);
		return c;
	}

	//////////////////////////////////////////////////////////////////////////
	// Boolean operators
	inline bool operator== (const mpcomplex& a, const mpcomplex& b) { return (mpc_cmp(a.mpc_srcptr(), b.mpc_srcptr()) ==0); }
	inline bool operator!= (const mpcomplex& a, const mpcomplex& b) { return (mpc_cmp(a.mpc_srcptr(), b.mpc_srcptr()) !=0); }

	inline bool operator== (const mpcomplex& a, const long int b) { return (mpc_cmp_si(a.mpc_srcptr(), b) ==0); }


	inline const mpreal real(const mpcomplex& x) { return x.real(); }
	inline const mpreal imag(const mpcomplex& x) { return x.imag(); }

	#define MPCOMPLEX_UNARY_MATH_FUNCTION_BODY(f)                 \
			mpcomplex y;          \
			mpc_##f(y.mpc_ptr(), x.mpc_srcptr(), MPC_RNDNN);           \
			return y; 

	inline const mpcomplex sqr(const mpcomplex& x) { MPCOMPLEX_UNARY_MATH_FUNCTION_BODY(sqr); }
	inline const mpcomplex sqrt(const mpcomplex& x) { MPCOMPLEX_UNARY_MATH_FUNCTION_BODY(sqrt); }

	inline const mpcomplex sin(const mpcomplex& x) { MPCOMPLEX_UNARY_MATH_FUNCTION_BODY(sin); }
	inline const mpcomplex cos(const mpcomplex& x) { MPCOMPLEX_UNARY_MATH_FUNCTION_BODY(cos); }
	inline const mpcomplex tan(const mpcomplex& x) { MPCOMPLEX_UNARY_MATH_FUNCTION_BODY(tan); }
	inline const mpcomplex cot(const mpcomplex& x) { return mpcomplex(1) / tan(x); }

	inline const mpcomplex sinh(const mpcomplex& x) { MPCOMPLEX_UNARY_MATH_FUNCTION_BODY(sinh); }
	inline const mpcomplex cosh(const mpcomplex& x) { MPCOMPLEX_UNARY_MATH_FUNCTION_BODY(cosh); }
	inline const mpcomplex tanh(const mpcomplex& x) { MPCOMPLEX_UNARY_MATH_FUNCTION_BODY(tanh); }
	inline const mpcomplex coth(const mpcomplex& x) { return mpcomplex(1) / tanh(x); }

	inline const mpcomplex asin(const mpcomplex& x) { MPCOMPLEX_UNARY_MATH_FUNCTION_BODY(asin); }
	inline const mpcomplex acos(const mpcomplex& x) { MPCOMPLEX_UNARY_MATH_FUNCTION_BODY(acos); }
	inline const mpcomplex atan(const mpcomplex& x) { MPCOMPLEX_UNARY_MATH_FUNCTION_BODY(atan); }

	inline const mpcomplex asinh(const mpcomplex& x) { MPCOMPLEX_UNARY_MATH_FUNCTION_BODY(asinh); }
	inline const mpcomplex acosh(const mpcomplex& x) { MPCOMPLEX_UNARY_MATH_FUNCTION_BODY(acosh); }
	inline const mpcomplex atanh(const mpcomplex& x) { MPCOMPLEX_UNARY_MATH_FUNCTION_BODY(atanh); }

	inline const mpcomplex conj(const mpcomplex& x) { MPCOMPLEX_UNARY_MATH_FUNCTION_BODY(conj); }
	inline const mpcomplex exp(const mpcomplex& x) { MPCOMPLEX_UNARY_MATH_FUNCTION_BODY(exp); }
	inline const mpcomplex log(const mpcomplex& x) { MPCOMPLEX_UNARY_MATH_FUNCTION_BODY(log); }

	inline const mpcomplex pow(const mpcomplex& x, long y)
	{
		mpcomplex r;
		mpc_pow_si(r.mpc_ptr(), x.mpc_srcptr(), y, MPC_RNDNN);
		return r;
	}

	inline const mpreal abs(const mpcomplex& x)
	{
		mpreal y;
		mpc_abs(y.mpfr_ptr(), x.mpc_srcptr(), mpreal::get_default_rnd());
		return y;
	}

	//////////////////////////////////////////////////////////////////////////
	// I/O
	inline std::ostream& mpcomplex::output(std::ostream& os) const 
	{
		std::ostringstream format;
		const std::ios::fmtflags flags = os.flags();

		format << ((flags & std::ios::showpos) ? "%+" : "%");
		if (os.precision() >= 0)
			format << '.' << os.precision() << "R*"
				   << ((flags & std::ios::floatfield) == std::ios::fixed ? 'f' :
					   (flags & std::ios::floatfield) == std::ios::scientific ? 'e' :
					   'g');
		else
			format << "R*e";

		char *s = NULL;
		if(!(mpfr_asprintf(&s, format.str().c_str(),
							mpfr::mpreal::get_default_rnd(),
							real().mpfr_srcptr())
			< 0))
		{
			os << std::string(s);
			mpfr_free_str(s);
		}
		
		if(!(mpfr_asprintf(&s, format.str().c_str(),
							mpfr::mpreal::get_default_rnd(),
							imag().mpfr_srcptr())
			< 0))
		{
			if (!(flags & std::ios::showpos) && (s[0] != '+') && (s[0] != '-'))
				os << " +" << std::string(s) << "i";
			else
				os << " " << std::string(s) << "i";
			mpfr_free_str(s);
		}
		return os;
	}

	inline std::ostream& operator<<(std::ostream& os, const mpcomplex& v)
	{
		return v.output(os);
	}

	inline std::string mpcomplex::toString(int n, int b) const
	{
		char* str = mpc_get_str(b, n, mpc_ptr(), MPC_RNDNN);
		std::string retVal(str);
		mpc_free_str(str);
		return retVal;
	}
}

#endif
