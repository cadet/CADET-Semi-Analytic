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

#include "DurbinsMethod.hpp"
#include "Constants.hpp"

namespace
{

template <typename func_t>
casema::util::SlicedVector<mpfr::mpreal> invertLaplaceImpl(const std::function<void(const mpfr::mpcomplex&, mpfr::mpcomplex*)>& f, int nOutputs, std::size_t nSummands, std::size_t precision, const mpfr::mpreal& tMax, const mpfr::mpreal& abscissa, mpfr::mpreal const* timeGrid, std::size_t nTime, func_t callback)
{
	casema::util::SlicedVector<mpfr::mpreal> result;
	result.reserve(nOutputs * nTime, nOutputs);
	for (int i = 0; i < nOutputs; ++i)
		result.pushBackSlice(nTime);

	const mpfr::mpreal T = tMax / 2 * mpfr::mpreal(1.01);

	std::vector<mpfr::mpreal> f0(nOutputs);

	{
		std::vector<mpfr::mpcomplex> temp(nOutputs);
		f(mpfr::mpcomplex(abscissa), temp.data());
		for (int i = 0; i < nOutputs; ++i)
			f0[i] = temp[i].real() / 2;
	}

	casema::util::SlicedVector<mpfr::mpcomplex> funcEvals;
	funcEvals.reserve(nOutputs * nSummands, nSummands);
	for (int i = 0; i < nSummands; ++i)
		funcEvals.pushBackSlice(nOutputs);

	std::size_t progress = 0;

	#pragma omp parallel
	{
#ifdef _OPENMP
		casema::Constants<mpfr::mpreal>::init(precision);
#endif

		#pragma omp for schedule(static)
		for (std::size_t k = 1; k <= nSummands; ++k)
		{
			const mpfr::mpreal kpit = k * casema::Constants<mpfr::mpreal>::pi() / T;
			f(mpfr::mpcomplex(abscissa, kpit), funcEvals[k-1]);

			#pragma omp critical
			{
				++progress;
				callback(progress / static_cast<double>(nSummands), 0);
			}
		}

		#pragma omp single nowait
		{
			progress = 0;
		}

		#pragma omp for schedule(static)
		for (std::size_t i = 0; i < nTime; ++i)
		{
			for (int j = 0; j < nOutputs; ++j)
				result(j, i) = f0[j];

			for (std::size_t k = 1; k <= nSummands; ++k)
			{
				const mpfr::mpreal kpit = k * casema::Constants<mpfr::mpreal>::pi() / T;

				mpfr::mpreal s;
				mpfr::mpreal c;
				mpfr::sin_cos(s, c, kpit * timeGrid[i]);
				for (int j = 0; j < nOutputs; ++j)
					result(j, i) += c * funcEvals(k-1, j).real() - s * funcEvals(k-1, j).imag();

				#pragma omp critical
				{
					++progress;
					callback(progress / static_cast<double>(nSummands * nTime), 1);
				}
			}

			for (int j = 0; j < nOutputs; ++j)
				result(j, i) *= exp(abscissa * timeGrid[i]) / T;
		}

#ifdef _OPENMP
		casema::Constants<mpfr::mpreal>::clear();
#endif
	}

	return result;
}

template <typename func_t>
casema::util::SlicedVector<mpfr::mpreal> invertLaplaceKahanImpl(const std::function<void(const mpfr::mpcomplex&, mpfr::mpcomplex*)>& f, int nOutputs, std::size_t nSummands, std::size_t precision, const mpfr::mpreal& tMax, const mpfr::mpreal& abscissa, mpfr::mpreal const* timeGrid, std::size_t nTime, func_t callback)
{
	casema::util::SlicedVector<mpfr::mpreal> result;
	result.reserve(nOutputs * nTime, nOutputs);
	for (int i = 0; i < nOutputs; ++i)
		result.pushBackSlice(nTime);

	const mpfr::mpreal T = tMax / 2 * mpfr::mpreal(1.01);

	std::vector<mpfr::mpreal> f0(nOutputs);

	{
		std::vector<mpfr::mpcomplex> temp(nOutputs);
		f(mpfr::mpcomplex(abscissa), temp.data());
		for (int i = 0; i < nOutputs; ++i)
			f0[i] = temp[i].real() / 2;
	}

	casema::util::SlicedVector<mpfr::mpcomplex> funcEvals;
	funcEvals.reserve(nOutputs * nSummands, nSummands);
	for (int i = 0; i < nSummands; ++i)
		funcEvals.pushBackSlice(nOutputs);

	std::size_t progress = 0;

	#pragma omp parallel
	{
#ifdef _OPENMP
		casema::Constants<mpfr::mpreal>::init(precision);
#endif

		#pragma omp for schedule(static)
		for (std::size_t k = 1; k <= nSummands; ++k)
		{
			const mpfr::mpreal kpit = k * casema::Constants<mpfr::mpreal>::pi() / T;
			f(mpfr::mpcomplex(abscissa, kpit), funcEvals[k-1]);

			#pragma omp critical
			{
				++progress;
				callback(progress / static_cast<double>(nSummands), 0);
			}
		}

		#pragma omp single nowait
		{
			progress = 0;
		}

		std::vector<mpfr::mpreal> runningError(nOutputs, 0.0);
		std::vector<mpfr::mpreal> temp(nOutputs, 0.0);
		std::vector<mpfr::mpreal> diff(nOutputs, 0.0);

		#pragma omp for schedule(static)
		for (std::size_t i = 0; i < nTime; ++i)
		{
			for (int j = 0; j < nOutputs; ++j)
			{
				runningError[j] = 0.0;
				result(j, i) = f0[j];
			}

			for (std::size_t k = 1; k <= nSummands; ++k)
			{
				const mpfr::mpreal kpit = k * casema::Constants<mpfr::mpreal>::pi() / T;

				mpfr::mpreal s;
				mpfr::mpreal c;
				mpfr::sin_cos(s, c, kpit * timeGrid[i]);

				// Kahan summation
				for (int j = 0; j < nOutputs; ++j)
				{
					diff[j] = c * funcEvals(k-1, j).real() - s * funcEvals(k-1, j).imag();
					diff[j] -= runningError[j];
					temp[j] = result(j, i);
					temp[j] += diff[j];
					runningError[j] = temp[j];
					runningError[j] -= result(j, i);
					runningError[j] -= diff[j];
					result(j, i) = temp[j];
				}

				#pragma omp critical
				{
					++progress;
					callback(progress / static_cast<double>(nSummands * nTime), 1);
				}
			}

			for (int j = 0; j < nOutputs; ++j)
				result(j, i) *= exp(abscissa * timeGrid[i]) / T;
		}

#ifdef _OPENMP
		casema::Constants<mpfr::mpreal>::clear();
#endif
	}

	return result;
}

void doNothing(double, int) { }

}

namespace casema
{
	util::SlicedVector<mpfr::mpreal> invertLaplace(const std::function<void(const mpfr::mpcomplex&, mpfr::mpcomplex*)>& f, int nOutputs, std::size_t nSummands, std::size_t precision, const mpfr::mpreal& tMax, const mpfr::mpreal& abscissa, mpfr::mpreal const* timeGrid, std::size_t nTime)
	{
		return invertLaplaceImpl(f, nOutputs, nSummands, precision, tMax, abscissa, timeGrid, nTime, doNothing);
	}

	util::SlicedVector<mpfr::mpreal> invertLaplaceKahan(const std::function<void(const mpfr::mpcomplex&, mpfr::mpcomplex*)>& f, int nOutputs, std::size_t nSummands, std::size_t precision, const mpfr::mpreal& tMax, const mpfr::mpreal& abscissa, mpfr::mpreal const* timeGrid, std::size_t nTime)
	{
		return invertLaplaceKahanImpl(f, nOutputs, nSummands, precision, tMax, abscissa, timeGrid, nTime, doNothing);
	}

	util::SlicedVector<mpfr::mpreal> invertLaplace(const std::function<void(const mpfr::mpcomplex&, mpfr::mpcomplex*)>& f, int nOutputs, std::size_t nSummands, std::size_t precision, const mpfr::mpreal& tMax, const mpfr::mpreal& abscissa, mpfr::mpreal const* timeGrid, std::size_t nTime, const std::function<void(double, int)>& callBack)
	{
		return invertLaplaceImpl(f, nOutputs, nSummands, precision, tMax, abscissa, timeGrid, nTime, callBack);
	}

	util::SlicedVector<mpfr::mpreal> invertLaplaceKahan(const std::function<void(const mpfr::mpcomplex&, mpfr::mpcomplex*)>& f, int nOutputs, std::size_t nSummands, std::size_t precision, const mpfr::mpreal& tMax, const mpfr::mpreal& abscissa, mpfr::mpreal const* timeGrid, std::size_t nTime, const std::function<void(double, int)>& callBack)
	{
		return invertLaplaceKahanImpl(f, nOutputs, nSummands, precision, tMax, abscissa, timeGrid, nTime, callBack);
	}
}
