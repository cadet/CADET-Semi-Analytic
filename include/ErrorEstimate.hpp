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

#ifndef CASEMA_ERROREST_HPP_
#define CASEMA_ERROREST_HPP_

#include <limits>
#include "ExpInt.hpp"
#include "Constants.hpp"

namespace casema
{

	namespace detail
	{

		template <typename real_t>
		real_t maxSplinePiece(const real_t& t0, const real_t& t1, const real_t& cons, const real_t& lin, const real_t& quad, const real_t& cub)
		{
			// Assume that spline piece is nonnegative, compute zeros of first derivative
			// We have f(t) = cub * (t-t0)^3 + quad * (t-t0)^2 + lin * (t-t0) + cons
			// Therefore, f'(t) = 3 * cub * (t-t0)^2 + 2 * quad * (t-t0) + lin
			// f'(t) = 0  <=>  t-t0 = (-quad +/- sqrt( quad^2 - 3 * lin * cub) ) / (3 * cub)
			
			// Stable solution
			const real_t radical = sqr(quad) - 3 * lin * cub;
			real_t maxVal = 0;
			if (radical > 0)
			{
				// Two real solutions
				
				real_t delta1;
				if (quad < 0)
					delta1 = (-quad + sqrt(radical)) / (3 * cub);
				else
					delta1 = (-quad - sqrt(radical)) / (3 * cub);

				const real_t delta2 = lin / (3 * cub * delta1);

				const real_t v1 = cons + delta1 * (lin + delta1 * (quad + delta1 * cub));
				const real_t v2 = cons + delta2 * (lin + delta2 * (quad + delta2 * cub));
				maxVal = max(real_t(0), max(v1, v2));
			}
			else if (abs(radical) <= std::numeric_limits<real_t>::epsilon())
			{
				// One real solution
				const real_t delta = -quad / (3 * cub);
				maxVal = max(real_t(0), cons + delta * (lin + delta * (quad + delta * cub)));
			}

			// Check boundaries
			const real_t delta = t1 - t0;
			const real_t right = cons + delta * (lin + delta * (quad + delta * cub));
			return max(maxVal, max(cons, right));
		}

		template <typename real_t, template <class T> class ModelType>
		real_t maxInletSpline(const ModelType<real_t>& model)
		{
			// Find maximum of nonnegative cubic spline
			real_t maxVal = 0;
			for (std::size_t i = 0; i < model.nInletSections; ++i)
			{
				maxVal = max(maxVal, maxSplinePiece(model.sectionTimes[i], model.sectionTimes[i+1], model.constCoeff[i], model.linCoeff[i], model.quadCoeff[i], model.cubicCoeff[i]));
			}
			return maxVal;
		}

		template <typename real_t>
		real_t inverseBoundedTruncError(const real_t& x, const real_t& tol = std::numeric_limits<real_t>::epsilon(), const unsigned int maxIter = 10000)
		{
			// Use Newton method

			// Construct initial guess
			const real_t one(1);
			real_t p = std::max(one, -log(x));

			// Calculate residual
			real_t residual = exp(-p) * (p + one) - x;

			// Perform Newton steps until accuracy is reached
			for (unsigned int i = 0; (i < maxIter) && ((std::abs(residual) > tol) || (std::abs(residual) > tol * std::abs(x))); ++i)
			{
				const real_t step = exp(p) / p * residual;

				// Stop if adding delta does not yield any change in p
				if (std::abs(step) <= machineEpsilon(p))
					break;

				// Try line search with at most 10 steps to decrease residual
				real_t alpha(1);
				real_t cand = p + alpha * step;
				real_t newResidual = exp(-cand) * (cand + one) - x;
				for (unsigned int i = 0; (i < 10) && (std::abs(newResidual) > std::abs(residual)); ++i)
				{
					alpha *= real_t(0.5);
					cand = p + alpha * step;
					newResidual = exp(-cand) * (cand + one) - x;
				}

				// Accept step
				p = cand;
				residual = newResidual;
			}

			return p;
		}

		template <typename real_t>
		real_t mInSplinePiece(const real_t& tU, const real_t& tL, const real_t& a, const real_t& cons, const real_t& lin, const real_t& quad, const real_t& cub)
		{
			const real_t six(6);
			const real_t two(2);
			const real_t dt = tU - tL;
			const real_t funcEval = abs(cons + dt * (lin + dt * (quad + dt * cub)));
			const real_t diffEval = abs(lin + dt * (two * quad + real_t(3) * dt * cub));
			const real_t diff2Eval = abs(two*quad + dt * six * cub);
			const real_t upperRest = ((six * abs(cub) / a + diff2Eval) / a + diffEval) / a;
			const real_t lowerRest = abs(cons) + (abs(lin) + ( two * abs(quad) + six * abs(cub) / a) / a) / a;

			return exp(-a * tU) * (funcEval + upperRest) + exp(-a * tL) * lowerRest;;
		}

		template <typename real_t, template <class T> class ModelType>
		real_t mInBoundStepLike(const real_t& a, const ModelType<real_t>& model)
		{
			real_t mIn(0);
			for (std::size_t i = 0; i < model.nInletSections; ++i)
			{
				mIn += mInSplinePiece(model.sectionTimes[i+1], model.sectionTimes[i], a, model.constCoeff[i], model.linCoeff[i], model.quadCoeff[i], model.cubicCoeff[i]);
			}
			return mIn;
		}


	}


	template <typename real_t, template <class T> class ModelType>
	real_t consistencyError(const real_t& a, const real_t& T, const ModelType<real_t>& model)
	{
		return detail::maxInletSpline(model) / expm1(a * T);
	}


	template <typename real_t, template <class T> class ModelType>
	real_t consistencyError(const real_t& a, const ModelType<real_t>& model)
	{
		return consistencyError(a, maxSimulationTime(model) * real_t(1.01), model);
	}


	template <typename real_t, template <class T> class ModelType, typename sum_t>
	real_t truncationErrorStep(const real_t& a, sum_t N, const real_t& T, const ModelType<real_t>& model)
	{
		const real_t Min = detail::mInBoundStepLike(a, model);
		const real_t sigma = sqrt(real_t(32)) * Min / Constants<real_t>::pi() * exp(model.colLength * model.velocity / (2 * model.colDispersion));
		const real_t gamma = model.colLength * sqrt(Constants<real_t>::pi() / (model.colDispersion * T));
		return exp(a * T) * sigma * expInt(gamma * sqrt(real_t(N+1)));
	}


	template <typename real_t, template <class T> class ModelType, typename sum_t>
	real_t truncationErrorBounded(const real_t& a, sum_t N, const real_t& T, const ModelType<real_t>& model)
	{
		const real_t Min = 2 * detail::maxInletSpline(model);
		const real_t sigma = sqrt(real_t(32)) * Min * exp(model.colLength * model.velocity / (2 * model.colDispersion));
		const real_t gamma = model.colLength * sqrt(Constants<real_t>::pi() / (model.colDispersion * T));
		const real_t sqrtNp1Gamma = sqrt(real_t(N+1)) * gamma;
		return exp(a * T) * sigma * exp(-sqrtNp1Gamma) * (sqrtNp1Gamma + 1) / sqr(gamma);
	}


	template <typename real_t, template <class T> class ModelType, typename sum_t>
	real_t truncationErrorStep(const real_t& a, sum_t N, const ModelType<real_t>& model)
	{
		return truncationErrorStep(a, N, maxSimulationTime(model) * real_t(1.01), model);
	}


	template <typename real_t, template <class T> class ModelType, typename sum_t>
	real_t truncationErrorBounded(const real_t& a, sum_t N, const ModelType<real_t>& model)
	{
		return truncationErrorBounded(a, N, maxSimulationTime(model) * real_t(1.01), model);
	}


	template <typename real_t, template <class T> class ModelType, typename sum_t>
	real_t truncationError(const real_t& a, sum_t N, const real_t& T, const ModelType<real_t>& model)
	{
		if (inletIsStep(model))
			return truncationErrorStep(a, N, T, model);
		else
			return truncationErrorBounded(a, N, T, model);
	}


	template <typename real_t, template <class T> class ModelType, typename sum_t>
	real_t truncationError(const real_t& a, sum_t N, const ModelType<real_t>& model)
	{
		return truncationError(a, N, maxSimulationTime(model) * real_t(1.01), model);
	}


	template <typename real_t, template <class T> class ModelType, typename sum_t>
	void paramsFromErrorStep(const real_t& abscissaSafety, const real_t& T, const ModelType<real_t>& model, const real_t& weight, const real_t& error, real_t& a, sum_t& N)
	{
		const real_t maxIn = detail::maxInletSpline(model);

		real_t nu = weight;
		if ((weight <= 0) || (weight >= 1))
			nu = real_t(1) / real_t(2);

		a = abscissaSafety + log(maxIn / (nu * error) + 1) / T;

		const real_t Min = detail::mInBoundStepLike(a, model);
		const real_t sigma = sqrt(real_t(32)) * Min / Constants<real_t>::pi() * exp(model.colLength * model.velocity / (2 * model.colDispersion));
		const real_t gamma = model.colLength * sqrt(Constants<real_t>::pi() / (model.colDispersion * T));

		N = static_cast<sum_t>(ceil(sqr(inverseExpInt((1-nu) * error / (sigma * exp(a * T))) / gamma)));
	}


	template <typename real_t, template <class T> class ModelType, typename sum_t>
	void paramsFromErrorStep(const real_t& abscissaSafety, const ModelType<real_t>& model, const real_t& weight, const real_t& error, real_t& a, sum_t& N)
	{
		paramsFromErrorStep(abscissaSafety, maxSimulationTime(model) * real_t(1.01), model, weight, error, a, N);
	}


	template <typename real_t, template <class T> class ModelType, typename sum_t>
	void paramsFromErrorBounded(const real_t& abscissaSafety, const real_t& T, const ModelType<real_t>& model, const real_t& weight, const real_t& error, real_t& a, sum_t& N)
	{
		const real_t maxIn = detail::maxInletSpline(model);
		const real_t Min = 2 * maxIn;
		const real_t sigma = sqrt(real_t(32)) * Min / Constants<real_t>::pi() * exp(model.colLength * model.velocity / (2 * model.colDispersion));
		const real_t gamma = model.colLength * sqrt(Constants<real_t>::pi() / (model.colDispersion * T));

		real_t nu = weight;
		if ((weight <= 0) || (weight >= 1))
			nu = real_t(1) / real_t(2);

		a = abscissaSafety + log(maxIn / (nu * error) + 1) / T;
		N = static_cast<sum_t>(ceil(sqr(detail::inverseBoundedTruncError((1-nu) * error * sqr(gamma) / (sigma * exp(a * T))) / gamma)));
	}


	template <typename real_t, template <class T> class ModelType, typename sum_t>
	void paramsFromErrorBounded(const real_t& abscissaSafety, const ModelType<real_t>& model, const real_t& weight, const real_t& error, real_t& a, sum_t& N)
	{
		paramsFromErrorBounded(abscissaSafety, maxSimulationTime(model) * real_t(1.01), model, weight, error, a, N);
	}


	template <typename real_t, template <class T> class ModelType, typename sum_t>
	void paramsFromError(const real_t& abscissaSafety, const real_t& T, const ModelType<real_t>& model, const real_t& weight, const real_t& error, real_t& a, sum_t& N)
	{
		if (inletIsStep(model))
			paramsFromErrorStep(abscissaSafety, T, model, weight, error, a, N);
		else
			paramsFromErrorBounded(abscissaSafety, T, model, weight, error, a, N);
	}


	template <typename real_t, template <class T> class ModelType, typename sum_t>
	void paramsFromError(const real_t& abscissaSafety, const ModelType<real_t>& model, const real_t& weight, const real_t& error, real_t& a, sum_t& N)
	{
		paramsFromError(abscissaSafety, maxSimulationTime(model) * real_t(1.01), model, weight, error, a, N);
	}
}

#endif
