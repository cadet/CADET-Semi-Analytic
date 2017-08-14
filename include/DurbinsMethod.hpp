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

#ifndef CASEMA_DURBIN_HPP_
#define CASEMA_DURBIN_HPP_

#include <vector>
#include <iostream>

#include "Constants.hpp"

namespace casema
{

	template <typename real_t, typename complex_t>
	class DurbinsMethod
	{
	public:

		DurbinsMethod() { }

		template <class lapFunc_t>
		std::vector<real_t> invert(const std::size_t precision, const std::size_t summands, const real_t& tMax, const real_t& sigma, const std::vector<real_t>& timeGrid, const std::size_t timeOffset, const lapFunc_t& f) const
		{
			const real_t T = tMax / real_t(2) * real_t(1.01);
			
			const real_t f0 = f(complex_t(sigma)).real() / real_t(2);
			std::vector<complex_t> funcEvals(summands);
			std::vector<real_t> result(timeGrid.size() - timeOffset);

			#pragma omp parallel
			{
				Constants<real_t>::init(precision);
				const real_t zero(0);

				#pragma omp single nowait
				{
					std::cout << "Evaluate Laplace solution at " << summands << " points" << std::endl;
				}

				#pragma omp for schedule(static)
				for (std::size_t k = 1; k <= summands; ++k)
				{
					const real_t kpit = k * Constants<real_t>::pi() / T;
					funcEvals[k-1] = f(sigma + complex_t(zero, kpit));
				}

				#pragma omp single nowait
				{
					std::cout << "Evaluate Fourier series for " << (timeGrid.size() - timeOffset) << " time points" << std::endl;
				}

				#pragma omp for schedule(static)
				for (std::size_t i = 0; i < timeGrid.size()-timeOffset; ++i)
				{
					const std::size_t timeIdx = i+timeOffset;

#ifdef CASEMA_KAHAN_SUMMATION
					real_t runningError = 0;
					real_t temp;
					real_t diff;
#endif

					result[i] = f0;
					for (std::size_t k = 1; k <= summands; ++k)
					{
						const real_t kpit = k * Constants<real_t>::pi() / T;

#ifdef CASEMA_KAHAN_SUMMATION
						// Kahan summation
						diff = cos(kpit * timeGrid[timeIdx]) * funcEvals[k-1].real() - sin(kpit * timeGrid[timeIdx]) * funcEvals[k-1].imag();
						diff -= runningError;
						temp = result[i];
						temp += diff;
						runningError = temp;
						runningError -= result[i];
						runningError -= diff;
						result[i] = std::move(temp); 
#else
						result[i] += cos(kpit * timeGrid[timeIdx]) * funcEvals[k-1].real() - sin(kpit * timeGrid[timeIdx]) * funcEvals[k-1].imag();
#endif
					}
					result[i] *= exp(sigma * timeGrid[timeIdx]) / T;
				}

				Constants<real_t>::clear();
			}
			return result;
		}

	protected:

	};

	template <typename real_t, typename complex_t>
	class PrecomputedDurbinsMethod
	{
	public:

		template <class lapFunc_t>
		PrecomputedDurbinsMethod(const std::size_t precision, const std::size_t summands, const real_t& tMax, const real_t& sigma, const lapFunc_t& f)
			: _precision(precision), _funcEvals(summands), _T(tMax / real_t(2) * real_t(1.01)), _sigma(sigma), _f0(f(complex_t(_sigma)).real() / real_t(2))
		{
			#pragma omp parallel
			{
				Constants<real_t>::init(precision);

				#pragma omp for schedule(static)
				for (std::size_t k = 1; k <= summands; ++k)
				{
					const real_t kpit = k * Constants<real_t>::pi() / _T;
					_funcEvals[k-1] = f(_sigma + complex_t(real_t(0), kpit));
				}

				Constants<real_t>::clear();
			}			
		}

		std::vector<real_t> invert(const std::vector<real_t>& timeGrid, const std::size_t timeOffset) const
		{
			std::vector<real_t> result(timeGrid.size() - timeOffset);

			#pragma omp parallel
			{
				Constants<real_t>::init(_precision);

				#pragma omp for schedule(static)
				for (std::size_t i = 0; i < timeGrid.size()-timeOffset; ++i)
				{
					const std::size_t timeIdx = i+timeOffset;

#ifdef CASEMA_KAHAN_SUMMATION
					real_t runningError = 0;
					real_t temp;
					real_t diff;
#endif

					result[i] = _f0;
					for (std::size_t k = 1; k <= _funcEvals.size(); ++k)
					{
						const real_t kpit = k * Constants<real_t>::pi() / _T;

#ifdef CASEMA_KAHAN_SUMMATION
						// Kahan summation
						diff = cos(kpit * timeGrid[timeIdx]) * _funcEvals[k-1].real() - sin(kpit * timeGrid[timeIdx]) * _funcEvals[k-1].imag();
						diff -= runningError;
						temp = result[i];
						temp += diff;
						runningError = temp;
						runningError -= result[i];
						runningError -= diff;
						result[i] = std::move(temp); 
#else
						result[i] += cos(kpit * timeGrid[timeIdx]) * _funcEvals[k-1].real() - sin(kpit * timeGrid[timeIdx]) * _funcEvals[k-1].imag();
#endif
					}
					result[i] *= exp(_sigma * timeGrid[timeIdx]) / _T;
				}

				Constants<real_t>::clear();
			}
			return result;
		}

		real_t operator()(const real_t& t) const { return invert(t); }

		real_t invert(const real_t& t) const
		{
#ifdef CASEMA_KAHAN_SUMMATION
			real_t runningError = 0;
			real_t temp;
			real_t diff;
#endif

			real_t result = _f0;
			for (std::size_t k = 1; k <= _funcEvals.size(); ++k)
			{
				const real_t kpit = k * Constants<real_t>::pi() / _T;

#ifdef CASEMA_KAHAN_SUMMATION
				// Kahan summation
				diff = cos(kpit * t) * _funcEvals[k-1].real() - sin(kpit * t) * _funcEvals[k-1].imag();
				diff -= runningError;
				temp = result;
				temp += diff;
				runningError = temp;
				runningError -= result;
				runningError -= diff;
				result = std::move(temp); 
#else
				result += cos(kpit * t) * _funcEvals[k-1].real() - sin(kpit * t) * _funcEvals[k-1].imag();
#endif
			}
			result *= exp(_sigma * t) / _T;
			return result;
		}

	protected:
		std::size_t _precision;
		std::vector<complex_t> _funcEvals;
		const real_t _T;
		const real_t _sigma;
		const real_t _f0;
	};

}

#endif
