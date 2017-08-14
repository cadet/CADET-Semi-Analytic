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

#ifndef CASEMA_EXTRAPOLATED_DURBIN_HPP_
#define CASEMA_EXTRAPOLATED_DURBIN_HPP_

#include <vector>
#include <iostream>

#include "Constants.hpp"
#include "CliExtrapolationParser.hpp"
#include "ExtrapolationHelpers.hpp"

namespace casema
{
	namespace detail
	{
		template <typename real_t>
		class VoidHookingPolicy
		{
		protected:
			void loopHook(const std::size_t summand, const real_t& seriesVal, const real_t& curEst, const real_t& a, const real_t& b) const { }
			void loopHook(const std::size_t summand, const real_t& seriesVal, const real_t& curEst, const real_t& a, const real_t& b, const real_t& c, const real_t& d, bool e, bool f) const { }
			void loopStartHook(bool combined) const { }
		};
	}

//	template <typename real_t, typename complex_t, template <typename T> class ExtrapolationPolicy, template <typename T> class HookingPolicy = detail::VoidHookingPolicy>
//	class ExtrapolatedDurbinsMethod : private ExtrapolationPolicy<real_t>, private HookingPolicy<real_t>
	template <typename real_t, typename complex_t, class ExtrapolationPolicy, template <typename T> class HookingPolicy = detail::VoidHookingPolicy>
	class ExtrapolatedDurbinsMethod : public ExtrapolationPolicy, public HookingPolicy<real_t>
	{
	public:

		ExtrapolatedDurbinsMethod() { }
		ExtrapolatedDurbinsMethod(const std::size_t precision, const std::size_t summands, const std::size_t minSummands, const real_t& tMax, const real_t& sigma)
			: _precision(precision), _summands(summands), _minSummands(minSummands), _tMax(tMax), _sigma(sigma) { }

		template <class lapFunc_t>
		std::vector<real_t> invert(const std::vector<real_t>& timeGrid, const std::size_t timeOffset, const lapFunc_t& f, real_t* const radius, std::size_t* const iterations) const
		{
			return invertCombined(timeGrid, timeOffset, f, radius, iterations);
		}

		template <class lapFunc_t>
		std::vector<real_t> invertSingle(const std::vector<real_t>& timeGrid, const std::size_t timeOffset, const lapFunc_t& f, real_t* const radius, std::size_t* const iterations) const
		{
			const real_t T = _tMax / real_t(2) * real_t(1.01);
			
			const real_t f0 = f(complex_t(_sigma)).real() / real_t(2);
			std::vector<complex_t> funcEvals(_summands);
			std::vector<real_t> result(timeGrid.size() - timeOffset);

			#pragma omp parallel
			{
				Constants<real_t>::init(_precision);
				auto sinExtrapolator = ExtrapolationPolicy::initializeThreaded(_summands);
				auto cosExtrapolator = ExtrapolationPolicy::initializeThreaded(_summands);

				const real_t zero(0);
				const real_t pit = Constants<real_t>::pi() / T;

				#pragma omp single nowait
				{
					std::cout << "Evaluate Laplace solution at " << _summands << " points" << std::endl;
				}

				#pragma omp for schedule(static)
				for (std::size_t k = 1; k <= _summands; ++k)
				{
					const real_t kpit = k * Constants<real_t>::pi() / T;
					funcEvals[k-1] = f(_sigma + complex_t(zero, kpit));
				}

				#pragma omp single nowait
				{
					std::cout << "Evaluate Fourier series for " << (timeGrid.size() - timeOffset) << " time points" << std::endl;
				}

				#pragma omp for schedule(static)
				for (std::size_t i = 0; i < timeGrid.size()-timeOffset; ++i)
				{
					const std::size_t timeIdx = i+timeOffset;
					result[i] = f0;

					// Extrapolate cos terms
					HookingPolicy<real_t>::loopStartHook(false);

					real_t cosSeriesVal = 0;
					real_t sinSeriesVal = 0;
					real_t cosLimit = 0;
					real_t sinLimit = 0;
					sinExtrapolator->reset();
					cosExtrapolator->reset();

#ifdef CASEMA_KAHAN_SUMMATION
					real_t cosRunningError = 0;
					real_t cosTemp;
					real_t cosDiff;

					real_t sinRunningError = 0;
					real_t sinTemp;
					real_t sinDiff;
#endif

					const real_t factor = exp(_sigma * timeGrid[timeIdx]) / T;

					for (std::size_t k = 1; k <= _summands; ++k)
					{
						const real_t kpit = k * pit;

						const bool cosConv = ExtrapolationPolicy::converged(cosExtrapolator);
						const bool sinConv = ExtrapolationPolicy::converged(sinExtrapolator);

						if (!cosConv || (k <= _minSummands))
						{
#ifdef CASEMA_KAHAN_SUMMATION
							// Kahan summation
							cosDiff = cos(kpit * timeGrid[timeIdx]) * funcEvals[k-1].real();
							cosDiff -= cosRunningError;
							cosTemp = cosSeriesVal;
							cosTemp += cosDiff;
							cosRunningError = cosTemp;
							cosRunningError -= cosSeriesVal;
							cosRunningError -= cosDiff;
							cosSeriesVal = std::move(cosTemp); 
#else
							cosSeriesVal += cos(kpit * timeGrid[timeIdx]) * funcEvals[k-1].real();
#endif
							cosLimit = cosExtrapolator->nextPoint(cosSeriesVal, real_t(1) / real_t(k));
						}
						if (!sinConv || (k <= _minSummands))
						{
#ifdef CASEMA_KAHAN_SUMMATION
							// Kahan summation
							sinDiff = sin(kpit * timeGrid[timeIdx]) * funcEvals[k-1].imag();
							sinDiff -= sinRunningError;
							sinTemp = sinSeriesVal;
							sinTemp += sinDiff;
							sinRunningError = sinTemp;
							sinRunningError -= sinSeriesVal;
							sinRunningError -= sinDiff;
							sinSeriesVal = std::move(sinTemp); 
#else
							sinSeriesVal += sin(kpit * timeGrid[timeIdx]) * funcEvals[k-1].imag();
#endif
							sinLimit = sinExtrapolator->nextPoint(sinSeriesVal, real_t(1) / real_t(k));
						}

						HookingPolicy<real_t>::loopHook(k, cosSeriesVal, cosLimit, sinSeriesVal, sinLimit, f0, factor, cosConv, sinConv);

						if (cosConv && sinConv && (k > _minSummands))
							break;
					}

					// Combine results to final value
					result[i] += cosLimit - sinLimit;
					result[i] *= factor;

					// Calculate precision estimate
					if (radius)
						radius[i] = factor * (ExtrapolationPolicy::radius(cosExtrapolator) + ExtrapolationPolicy::radius(sinExtrapolator));

					if (iterations)
						iterations[i] = std::max(sinExtrapolator->iterations(), cosExtrapolator->iterations());
				}

				ExtrapolationPolicy::cleanUpThreaded(cosExtrapolator);
				ExtrapolationPolicy::cleanUpThreaded(sinExtrapolator);
				Constants<real_t>::clear();
			}
			return result;
		}

		template <class lapFunc_t>
		std::vector<real_t> invertCombined(const std::vector<real_t>& timeGrid, const std::size_t timeOffset, const lapFunc_t& f, real_t* const radius, std::size_t* const iterations) const
		{
			const real_t T = _tMax / real_t(2) * real_t(1.01);
			
			const real_t f0 = f(complex_t(_sigma)).real() / real_t(2);
			std::vector<complex_t> funcEvals(_summands);
			std::vector<real_t> result(timeGrid.size() - timeOffset);

			#pragma omp parallel
			{
				Constants<real_t>::init(_precision);
				auto extrapolator = ExtrapolationPolicy::initializeThreaded(_summands);

				const real_t zero(0);
				const real_t pit = Constants<real_t>::pi() / T;

				#pragma omp single nowait
				{
					std::cout << "Evaluate Laplace solution at " << _summands << " points" << std::endl;
				}

				#pragma omp for schedule(static)
				for (std::size_t k = 1; k <= _summands; ++k)
				{
					const real_t kpit = k * Constants<real_t>::pi() / T;
					funcEvals[k-1] = f(_sigma + complex_t(zero, kpit));
				}

				#pragma omp single nowait
				{
					std::cout << "Evaluate Fourier series for " << (timeGrid.size() - timeOffset) << " time points" << std::endl;
				}

				#pragma omp for schedule(static)
				for (std::size_t i = 0; i < timeGrid.size()-timeOffset; ++i)
				{
					const std::size_t timeIdx = i+timeOffset;
					result[i] = f0;

					// Extrapolate
					HookingPolicy<real_t>::loopStartHook(true);

					real_t seriesVal = 0;
					real_t limit = 0;
					extrapolator->reset();
					const real_t factor = exp(_sigma * timeGrid[timeIdx]) / T;

#ifdef CASEMA_KAHAN_SUMMATION
					real_t runningError = 0;
					real_t temp;
					real_t diff;
#endif

					for (std::size_t k = 1; k <= _summands; ++k)
					{
						const real_t kpit = k * pit;

#ifdef CASEMA_KAHAN_SUMMATION
						// Kahan summation
						diff = cos(kpit * timeGrid[timeIdx]) * funcEvals[k-1].real() - sin(kpit * timeGrid[timeIdx]) * funcEvals[k-1].imag();
						diff -= runningError;
						temp = seriesVal;
						temp += diff;
						runningError = temp;
						runningError -= seriesVal;
						runningError -= diff;
						seriesVal = std::move(temp); 
#else
						seriesVal += cos(kpit * timeGrid[timeIdx]) * funcEvals[k-1].real() - sin(kpit * timeGrid[timeIdx]) * funcEvals[k-1].imag();
#endif
						limit = extrapolator->nextPoint(seriesVal, real_t(1) / real_t(k));

						HookingPolicy<real_t>::loopHook(k, seriesVal, limit, f0, factor);

						if (ExtrapolationPolicy::converged(extrapolator) && (k > _minSummands))
							break;
					}
					const real_t rad = ExtrapolationPolicy::radius(extrapolator);

					// Combine results to final value
					result[i] += limit;
					result[i] *= factor;

					if (iterations)
						iterations[i] = extrapolator->iterations();

					// Calculate precision estimate
					if (radius)
						radius[i] = factor * rad;
				}

				ExtrapolationPolicy::cleanUpThreaded(extrapolator);
				Constants<real_t>::clear();
			}
			return result;
		}

	protected:
		std::size_t _precision;
		std::size_t _summands;
		std::size_t _minSummands;
		real_t _tMax;
		real_t _sigma;
	};


	template <typename real_t>
	class ConsensusExtrapolatorPolicy
	{
	public:

		struct ExtrapolatorOptions
		{
			std::size_t convTimes;
			real_t convThreshold;
			std::size_t consAgree;
			std::string extraMethods;
		};

		template <typename Options_t>
		static ExtrapolatorOptions convertOptions(const Options_t& opts)
		{
			ExtrapolatorOptions convOpts;

			convOpts.convTimes = opts.convTimes;
			convOpts.convThreshold = opts.convThreshold;
			convOpts.consAgree = opts.consAgree;
			convOpts.extraMethods = opts.extraMethods;

			return convOpts;
		}

		ConsensusExtrapolatorPolicy() { }

		template <typename Options_t>
		ConsensusExtrapolatorPolicy(const Options_t& opts) : _opts(convertOptions(opts)) { }

		template <typename Options_t>
		void options(const Options_t& opts) { _opts = convertOptions(opts); }

		const ExtrapolatorOptions& options() const { return _opts; }

	protected:

		casema::ConsensusEstimator<real_t>* initializeThreaded(std::size_t summands) const
		{
			casema::ConsensusEstimator<real_t>* estimator = new casema::ConsensusEstimator<real_t>();
			casema::addExtrapolationMethods(_opts.extraMethods, *estimator, _opts);
			return estimator;
		}

		void cleanUpThreaded(casema::ConsensusEstimator<real_t>* const estimator) const
		{
			delete estimator;
		}

		const real_t radius(casema::ConsensusEstimator<real_t> const* const estimator) const { return estimator->radius(); }

		bool converged(casema::ConsensusEstimator<real_t> const* const estimator) const
		{
			return false;
		}

		ExtrapolatorOptions _opts;
	};


	template <typename real_t, template <class T> class Extrapolator_t>
	class SingleExtrapolatorPolicy
	{
	public:

		struct ExtrapolatorOptions
		{
			std::size_t convTimes;
			real_t convThresholdAbs;
			real_t convThresholdRel;
		};

		template <typename Options_t>
		static ExtrapolatorOptions convertOptions(const Options_t& opts)
		{
			ExtrapolatorOptions convOpts;

			convOpts.convTimes = opts.convTimes;
			convOpts.convThresholdAbs = opts.convThresholdAbs;
			convOpts.convThresholdRel = opts.convThresholdRel;

			return convOpts;
		}

		SingleExtrapolatorPolicy() { }

		template <typename Options_t>
		SingleExtrapolatorPolicy(const Options_t& opts) : _opts(convertOptions(opts)) { }

		template <typename Options_t>
		void options(const Options_t& opts) { _opts = convertOptions(opts); }

		const ExtrapolatorOptions& options() const { return _opts; }

	protected:

		casema::ConvergenceMonitor<real_t, Extrapolator_t>* initializeThreaded(std::size_t summands) const
		{
			casema::ConvergenceMonitor<real_t, Extrapolator_t>* estimator = new casema::ConvergenceMonitor<real_t, Extrapolator_t>();
			estimator->thresholdAbs(_opts.convThresholdAbs);
			estimator->thresholdRel(_opts.convThresholdRel);
			estimator->times(_opts.convTimes);
			return estimator;
		}

		void cleanUpThreaded(casema::ConvergenceMonitor<real_t, Extrapolator_t>* const estimator) const
		{
			delete estimator;
		}

		const real_t radius(casema::ConvergenceMonitor<real_t, Extrapolator_t> const* const estimator) const { return 0; }

		bool converged(casema::ConvergenceMonitor<real_t, Extrapolator_t> const* const estimator) const
		{
			return estimator->converged();
		}

		ExtrapolatorOptions _opts;
	};

/*
	template <template <class T> class Extrapolator_t>
	struct SingleExtrapolatorPolicyGenerator
	{
		template <typename real_t> using PolicyType = SingleExtrapolatorPolicy<real_t, Extrapolator_t>;
	};
*/

}

#endif
