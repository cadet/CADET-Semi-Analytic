
#ifndef CASEMA_MPAD_HPP_
#define CASEMA_MPAD_HPP_

#include <mpreal.h>
#include "MPComplex.hpp"

#ifdef CASEMA_USE_FADBAD

#include <fadbad.h>

namespace fadbad
{
	template <> struct Op<mpfr::mpreal>
	{
		typedef mpfr::mpreal Base;
		static Base myInteger(const int i) { return Base(i); }
		static Base myZero() { return Base(0); }
		static Base myOne() { return Base(1);}
		static Base myTwo() { return Base(2); }
		static Base myPI() { return mpfr::const_pi(); }
		static Base myPos(const Base& x) { return +x; }
		static Base myNeg(const Base& x) { return -x; }
		template <typename U> static Base& myCadd(Base& x, const U& y) { return x+=y; }
		template <typename U> static Base& myCsub(Base& x, const U& y) { return x-=y; }
		template <typename U> static Base& myCmul(Base& x, const U& y) { return x*=y; }
		template <typename U> static Base& myCdiv(Base& x, const U& y) { return x/=y; }
		static Base myInv(const Base& x) { return myOne()/x; }
		static Base mySqr(const Base& x) { return ::mpfr::sqr(x); }
		template <typename X, typename Y>
		static Base myPow(const X& x, const Y& y) { return ::mpfr::pow(x,y); }
		static Base mySqrt(const Base& x) { return ::mpfr::sqrt(x); }
		static Base myLog(const Base& x) { return ::mpfr::log(x); }
		static Base myExp(const Base& x) { return ::mpfr::exp(x); }
		static Base mySin(const Base& x) { return ::mpfr::sin(x); }
		static Base myCos(const Base& x) { return ::mpfr::cos(x); }
		static Base myTan(const Base& x) { return ::mpfr::tan(x); }
		static Base myCot(const Base& x) { return ::mpfr::cot(x); }
		static Base myAsin(const Base& x) { return ::mpfr::asin(x); }
		static Base myAcos(const Base& x) { return ::mpfr::acos(x); }
		static Base myAtan(const Base& x) { return ::mpfr::atan(x); }
		static Base mySinh(const Base& x) { return ::mpfr::sinh(x); }
		static Base myCosh(const Base& x) { return ::mpfr::cosh(x); }
		static Base myTanh(const Base& x) { return ::mpfr::tanh(x); }
		static Base myCoth(const Base& x) { return ::mpfr::coth(x); }
		static bool myEq(const Base& x, const Base& y) { return x==y; }
		static bool myNe(const Base& x, const Base& y) { return x!=y; }
		static bool myLt(const Base& x, const Base& y) { return x<y; }
		static bool myLe(const Base& x, const Base& y) { return x<=y; }
		static bool myGt(const Base& x, const Base& y) { return x>y; }
		static bool myGe(const Base& x, const Base& y) { return x>=y; }
	};

	template <> struct Op<mpfr::mpcomplex>
	{
		typedef mpfr::mpcomplex Base;
		static Base myInteger(const int i) { return Base(i); }
		static Base myZero() { return Base(0); }
		static Base myOne() { return Base(1);}
		static Base myTwo() { return Base(2); }
		static Base myPI() { return mpfr::const_pi(); }
		static Base myPos(const Base& x) { return x; }
		static Base myNeg(const Base& x) { return -x; }
		template <typename U> static Base& myCadd(Base& x, const U& y) { return x+=y; }
		template <typename U> static Base& myCsub(Base& x, const U& y) { return x-=y; }
		template <typename U> static Base& myCmul(Base& x, const U& y) { return x*=y; }
		template <typename U> static Base& myCdiv(Base& x, const U& y) { return x/=y; }
		static Base myInv(const Base& x) { return myOne()/x; }
		static Base mySqr(const Base& x) { return ::mpfr::sqr(x); }
		template <typename X, typename Y>
		static Base myPow(const X& x, const Y& y) { return ::mpfr::pow(x,y); }
		static Base mySqrt(const Base& x) { return ::mpfr::sqrt(x); }
		static Base myLog(const Base& x) { return ::mpfr::log(x); }
		static Base myExp(const Base& x) { return ::mpfr::exp(x); }
		static Base mySin(const Base& x) { return ::mpfr::sin(x); }
		static Base myCos(const Base& x) { return ::mpfr::cos(x); }
		static Base myTan(const Base& x) { return ::mpfr::tan(x); }
		static Base myCot(const Base& x) { return ::mpfr::cot(x); }
		static Base myAsin(const Base& x) { return ::mpfr::asin(x); }
		static Base myAcos(const Base& x) { return ::mpfr::acos(x); }
		static Base myAtan(const Base& x) { return ::mpfr::atan(x); }
		static Base mySinh(const Base& x) { return ::mpfr::sinh(x); }
		static Base myCosh(const Base& x) { return ::mpfr::cosh(x); }
		static Base myTanh(const Base& x) { return ::mpfr::tanh(x); }
		static Base myCoth(const Base& x) { return ::mpfr::coth(x); }
		static bool myEq(const Base& x, const Base& y) { return x==y; }
		static bool myNe(const Base& x, const Base& y) { return x!=y; }
		// Complex values do not possess an order
		static bool myLt(const Base& x, const Base& y) { return false; }
		static bool myLe(const Base& x, const Base& y) { return false; }
		static bool myGt(const Base& x, const Base& y) { return false; }
		static bool myGe(const Base& x, const Base& y) { return false; }
	};

}

#endif

#include <cppad/configure.hpp>
#include <cppad/base_require.hpp>

namespace CppAD 
{
	inline mpfr::mpreal CondExpOp( 
		enum CompareOp     cop          ,
		const mpfr::mpreal&       left         ,
		const mpfr::mpreal&       right        , 
		const mpfr::mpreal&       exp_if_true  , 
		const mpfr::mpreal&       exp_if_false )
	{	return CondExpTemplate(cop, left, right, exp_if_true, exp_if_false);
	}

	CPPAD_COND_EXP_REL(mpfr::mpreal)


	inline bool EqualOpSeq(const mpfr::mpreal& x, const mpfr::mpreal& y)
	{	return x == y; }


	inline bool IdenticalPar(const mpfr::mpreal& x)
	{	return true; }
	inline bool IdenticalZero(const mpfr::mpreal& x)
	{	return mpfr::iszero(x); }
	inline bool IdenticalOne(const mpfr::mpreal& x)
	{	return (x == mpfr::mpreal("1")); }
	inline bool IdenticalEqualPar(const mpfr::mpreal& x, const mpfr::mpreal& y)
	{	return (x == y); }


	inline int Integer(const mpfr::mpreal& x)
	{	return static_cast<int>(x.toLong()); }

	CPPAD_AZMUL( mpfr::mpreal )

	inline bool GreaterThanZero(const mpfr::mpreal& x)
	{	return x > 0.; }
	inline bool GreaterThanOrZero(const mpfr::mpreal& x)
	{	return x >= 0.; }
	inline bool LessThanZero(const mpfr::mpreal& x)
	{	return x < 0.; }
	inline bool LessThanOrZero(const mpfr::mpreal& x)
	{	return x <= 0.; }
	inline bool abs_geq(const mpfr::mpreal& x, const mpfr::mpreal& y)
	{	return mpfr::cmpabs(x, y) >= 0; }

#define CPPAD_MPREAL_MATH_UNARY(Type, Fun) \
	using ::mpfr::Fun;


	CPPAD_MPREAL_MATH_UNARY(mpfr::mpreal, acos)
	CPPAD_MPREAL_MATH_UNARY(mpfr::mpreal, asin)
	CPPAD_MPREAL_MATH_UNARY(mpfr::mpreal, atan)
	CPPAD_MPREAL_MATH_UNARY(mpfr::mpreal, cos)
	CPPAD_MPREAL_MATH_UNARY(mpfr::mpreal, cosh)
	CPPAD_MPREAL_MATH_UNARY(mpfr::mpreal, exp)
	CPPAD_MPREAL_MATH_UNARY(mpfr::mpreal, fabs)
	CPPAD_MPREAL_MATH_UNARY(mpfr::mpreal, log)
	CPPAD_MPREAL_MATH_UNARY(mpfr::mpreal, log10)
	CPPAD_MPREAL_MATH_UNARY(mpfr::mpreal, sin)
	CPPAD_MPREAL_MATH_UNARY(mpfr::mpreal, sinh)
	CPPAD_MPREAL_MATH_UNARY(mpfr::mpreal, sqrt)
	CPPAD_MPREAL_MATH_UNARY(mpfr::mpreal, tan)
	CPPAD_MPREAL_MATH_UNARY(mpfr::mpreal, tanh)
	CPPAD_MPREAL_MATH_UNARY(mpfr::mpreal, coth)

	inline mpfr::mpreal sign(const mpfr::mpreal& x)
	{	
		return (mpfr::iszero(x) ? 0 : mpfr::sgn(x));
	}


	template <>
	class numeric_limits<mpfr::mpreal> {
	public:
		// machine epsilon
		static mpfr::mpreal epsilon(void)
		{	return std::numeric_limits<mpfr::mpreal>::epsilon(); }
		// minimum positive normalized value
		static mpfr::mpreal min(void)
		{	return std::numeric_limits<mpfr::mpreal>::min(); }
		// maximum finite value
		static mpfr::mpreal max(void)
		{	return std::numeric_limits<mpfr::mpreal>::max(); }
		static mpfr::mpreal quiet_NaN(void)
		{	return std::numeric_limits<mpfr::mpreal>::quiet_NaN(); }
	};
	// deprecated machine epsilon
	template <> 
	inline mpfr::mpreal epsilon<mpfr::mpreal>(void)
	{	return numeric_limits<mpfr::mpreal>::epsilon(); }
}

#endif
