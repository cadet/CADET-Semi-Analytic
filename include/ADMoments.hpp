// =============================================================================
//  CADET-semi-analytic - The semi analytic extension of
//  		CADET - The Chromatography Analysis and Design Toolkit
//  
//  Copyright © 2015-2018: Samuel Leweke¹
//                                      
//    ¹ Forschungszentrum Juelich GmbH, IBG-1, Juelich, Germany.
//  
//  All rights reserved. This program and the accompanying materials
//  are made available under the terms of the GNU Public License v3.0 (or, at
//  your option, any later version) which accompanies this distribution, and
//  is available at http://www.gnu.org/licenses/gpl.html
// =============================================================================

#ifndef CASEMA_ADMOMENTS_HPP_
#define CASEMA_ADMOMENTS_HPP_

#include <vector>

#include "MPAD.hpp"

#ifdef CASEMA_USE_FADBAD
	#include <tadiff.h>
#endif

#include <cppad/cppad.hpp>

namespace casema
{

	template <class real_t>
	class MomentGenerator
	{
	public:
		MomentGenerator() { }
		virtual ~MomentGenerator() { }

		virtual real_t moment(const real_t& pos, std::size_t order) = 0;
		virtual std::vector<real_t> moments(const real_t& pos, std::size_t order) = 0;

		virtual real_t moment(const real_t& pos, const real_t& normalization, std::size_t order)
		{
			const real_t mom = moment(pos, order);

			if (order >= 1)
				return mom / normalization;

			return mom;
		}

		virtual std::vector<real_t> moments(const real_t& pos, const real_t& normalization, std::size_t order)
		{
			std::vector<real_t> mom = moments(pos, order);
			for (std::size_t i = 1; i < mom.size(); ++i)
				mom[i] /= normalization;
			return mom;
		}

	protected:

		template <typename T>
		real_t convertDerivativeToMoment(const T& dfdx, std::size_t order) const
		{
			real_t fact(1);
			for (std::size_t i = 2; i < order+1; ++i)
			{
				fact *= i;
			}

			if (order % 2 == 1)
				fact = -fact;

			return fact * dfdx[order];
		}

		template <typename T>
		std::vector<real_t> convertDerivativeToMoments(const T& dfdx, std::size_t order) const 
		{
			std::vector<real_t> output(order+1);
			real_t fact(1);
			for (std::size_t i = 0; i < order+1; ++i)
			{
				if (i > 1)
					fact *= i;

				if (i % 2 == 1)
					output[i] = -fact * dfdx[i];
				else
					output[i] = fact * dfdx[i];
			}
			return output;
		}
	};

#ifdef CASEMA_USE_FADBAD

	template <class real_t, class Solution_t>
	class FadBadMomentGenerator : public MomentGenerator<real_t>
	{
	public:
		FadBadMomentGenerator(Solution_t func) : _func(func) { }
		virtual ~FadBadMomentGenerator() { }

		virtual real_t moment(const real_t& pos, std::size_t order)
		{
			// Taylor expansion
			fadbad::T<real_t> x;
			x = pos;
			x[1] = real_t(1);

			fadbad::T<real_t> val = _func(x);
			val.eval(order);
			return this->convertDerivativeToMoment(val, order);
		}

		virtual std::vector<real_t> moments(const real_t& pos, std::size_t order)
		{
			// Taylor expansion
			fadbad::T<real_t> x;
			x = pos;
			x[1] = real_t(1);

			fadbad::T<real_t> val = _func(x);
			val.eval(order);
			return this->convertDerivativeToMoments(val, order);
		}

	protected:
		// TODO: Find out why FADBAD is producing an access violation
		//       if we have const Solution_t& here
		Solution_t _func;
	};

#endif

	template <class real_t, class Solution_t>
	class CppADMomentGenerator : public MomentGenerator<real_t>
	{
	public:
		CppADMomentGenerator(const Solution_t& func) : _adFunc(nullptr), _func(func) { setup(); }
		virtual ~CppADMomentGenerator() { if (_adFunc) delete _adFunc; }

		virtual real_t moment(const real_t& evalPoint, std::size_t order)
		{
			// Where to evaluate
			std::vector<real_t> pos(order+1, real_t(0));
			pos[0] = evalPoint;
			pos[1] = real_t(1);

			const std::vector<real_t> dfdx = _adFunc->Forward(order, pos);
			return this->convertDerivativeToMoment(dfdx, order);
		}

		virtual std::vector<real_t> moments(const real_t& evalPoint, std::size_t order)
		{
			// Where to evaluate
			std::vector<real_t> pos(order+1, real_t(0));
			pos[0] = evalPoint;
			pos[1] = real_t(1);

			const std::vector<real_t> dfdx = _adFunc->Forward(order, pos);
			return this->convertDerivativeToMoments(dfdx, order);
		}

	protected:

		void setup()
		{
			if (_adFunc)
				return;

			// Evaluation point does not matter since no branching is involved
			std::vector< CppAD::AD<real_t> > x(1);
			x[0] = real_t(1);

			// Start recording
			CppAD::Independent(x);
			// Execute function
			std::vector< CppAD::AD<real_t> > y(1);
			y[0] = _func(x[0]);
			// Stop recording
			_adFunc = new CppAD::ADFun<real_t>(x, y);
		}

		CppAD::ADFun<real_t>* _adFunc;
		const Solution_t& _func;
	};

}

#endif
