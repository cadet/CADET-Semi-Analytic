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

#include "model/GeneralRateModel2D.hpp"
#include "io/ParameterProvider.hpp"
#include "model/ParamReaderScopes.hpp"
#include "model/ParamReaderUtils.hpp"
#include "Exceptions.hpp"
#include "Constants.hpp"
#include "CompileTimeConfig.hpp"

#include <unordered_map>

namespace casema
{

namespace model
{

/**
 * @brief Creates a ConvectionDispersionLikeModel
 */
GeneralRateModel2D::GeneralRateModel2D(int unitOpIdx) : ColumnWithPoreDiffusion(unitOpIdx), _crossAreaZone(0), _j0nSq(0), _besselZerosPhi(0)
{
}

GeneralRateModel2D::~GeneralRateModel2D() CASEMA_NOEXCEPT
{
}

/**
 * @brief Reads model parameters
 * @details Only reads parameters that do not affect model structure (e.g., discretization).
 * @param [in] paramProvider Parameter provider for reading parameters
 * @return @c true if configuration went fine, @c false otherwise
 */
bool GeneralRateModel2D::configure(io::IParameterProvider& paramProvider)
{
	ColumnWithPoreDiffusion::configure(paramProvider);

	_colRadius = paramProvider.getDouble("COL_RADIUS");
	_crossSection = Constants<mpfr::mpreal>::pi() * sqr(_colRadius);

	paramProvider.pushScope("discretization");

	const int nRad = paramProvider.getInt("NRAD");

	const std::string rdt = paramProvider.getString("RADIAL_DISC_TYPE");
	if (rdt == "EQUIVOLUME")
	{
		_radialEdges.resize(nRad+1);
		const mpfr::mpreal volPerCompartment = _colRadius * _colRadius / nRad;
		const mpfr::mpreal pi = Constants<mpfr::mpreal>::pi();

		_radialEdges[0] = 0.0;
		for (unsigned int r = 0; r < nRad; ++r)
		{
			// Set last edge to _colRadius for exact geometry
			if (r == nRad - 1)
				_radialEdges[r+1] = _colRadius;
			else
				_radialEdges[r+1] = sqrt(volPerCompartment + _radialEdges[r] * _radialEdges[r]);
		}
	}
	else if (rdt == "USER_DEFINED")
	{
		const std::vector<double> radEdges = paramProvider.getDoubleArray("RADIAL_COMPARTMENTS");

		if (radEdges.size() < nRad + 1)
			throw InvalidParameterException("Number of elements in field RADIAL_COMPARTMENTS is less than NRAD + 1 (" + std::to_string(nRad + 1) + ")");

		_radialEdges = std::move(toMPreal(radEdges));
	}
	else
	{
		_radialEdges.resize(nRad+1);
		const mpfr::mpreal h = _colRadius / nRad;
		const mpfr::mpreal pi = Constants<mpfr::mpreal>::pi();

		_radialEdges[0] = 0.0;
		for (unsigned int r = 0; r < nRad; ++r)
		{
			// Set last edge to _colRadius for exact geometry
			if (r == nRad - 1)
				_radialEdges[r+1] = _colRadius;
			else
				_radialEdges[r+1] = h * (r + 1);
		}
	}

	paramProvider.popScope();

	const std::vector<double> colDispersionRadial = paramProvider.getDoubleArray("COL_DISPERSION_RADIAL");
	if ((colDispersionRadial.size() != 1) && (colDispersionRadial.size() != nRad))
		throw InvalidParameterException("Field COL_DISPERSION_RADIAL has to have 1 or " + std::to_string(nRad) + " elements");

	for (int i = 1; i < colDispersionRadial.size(); ++i)
	{
		if (colDispersionRadial[0] != colDispersionRadial[i])
			throw InvalidParameterException("Field COL_DISPERSION_RADIAL has to have the same value if it is an array");
	}

	_colDispersionRadial = colDispersionRadial[0];

	_crossAreaZone = std::vector<mpfr::mpreal>(nRad);
	for (int i = 0; i < nRad; ++i)
		_crossAreaZone[i] = 2 / (sqr(_radialEdges[i+1] / _colRadius) - sqr(_radialEdges[i] / _colRadius));

	return true;
}

void GeneralRateModel2D::setVelocityFromFlowRate(mpfr::mpreal const* in)
{
	const mpfr::mpreal pi = Constants<mpfr::mpreal>::pi();
	_velocity = in[0] / (_colPorosity * pi * (sqr(_radialEdges[1]) - sqr(_radialEdges[0])));

	for (int i = 1; i < _radialEdges.size() - 1; ++i)
	{
		const mpfr::mpreal v = in[i] / (_colPorosity * pi * (sqr(_radialEdges[i+1]) - sqr(_radialEdges[i])));
		if (abs(v - _velocity) >= 1e-14)
			throw InvalidParameterException("Flow rates in radial zones are not equal (zone " + std::to_string(i) + ", diff " + abs(v - _velocity).toString(6) + ")");
	}
}

void GeneralRateModel2D::setBesselZeros(mpfr::mpreal const* zeros, int n) CASEMA_NOEXCEPT
{
	const mpfr::mpreal factor = _colDispersionRadial / sqr(_colRadius);

	_besselZerosPhi.clear();
	_besselZerosPhi.reserve(n);
	for (int i = 0; i < n; ++i)
		_besselZerosPhi.push_back(sqr(zeros[i]) * factor);

	_j0nSq.clear();
	_j0nSq.reserve(n);
	_j0nSq.push_back(2);
	for (int i = 1; i < n; ++i)
		_j0nSq.push_back(2 / sqr(mpfr::besselj0(zeros[i])));

	const int nZones = _radialEdges.size() - 1;
	_besselCoeffs.clear();
	_besselCoeffs.reserve(n * nZones, nZones);

#ifdef ENABLE_BESSEL_TRUNCATION
	const mpfr::mpreal eps = 10 * std::numeric_limits<mpfr::mpreal>::epsilon();
#endif

	for (int i = 0; i < nZones; ++i)
	{
		const mpfr::mpreal& rho_mp1 = _radialEdges[i+1] / _colRadius;
		const mpfr::mpreal& rho_m = _radialEdges[i] / _colRadius;

		_besselCoeffs.pushBack((sqr(rho_mp1) - sqr(rho_m))/2);
		for (int j = 1; j < n; ++j)
		{
#ifdef ENABLE_BESSEL_TRUNCATION
			const mpfr::mpreal a = mpfr::besselj1(zeros[j] * rho_mp1);
			const mpfr::mpreal b = mpfr::besselj1(zeros[j] * rho_m);
			const mpfr::mpreal a2 = (abs(a) <= eps) ? 0 : a;
			const mpfr::mpreal b2 = (abs(b) <= eps) ? 0 : b;
			_besselCoeffs.pushBackInLastSlice((rho_mp1 * a2 - rho_m * b2) / zeros[j]);
#else
			_besselCoeffs.pushBackInLastSlice((rho_mp1 * mpfr::besselj1(zeros[j] * rho_mp1) - rho_m * mpfr::besselj1(zeros[j] * rho_m)) / zeros[j]);
#endif
		}
	}
}

mpfr::mpreal GeneralRateModel2D::timeDomainUpperBound(mpfr::mpreal const* in) const CASEMA_NOEXCEPT
{
	mpfr::mpreal outMax(0);
	for (int i = 0; i < _radialEdges.size()-1; ++i)
		outMax = max(outMax, in[i]);

	return outMax + _initC;
}

void GeneralRateModel2D::evaluate(const mpfr::mpcomplex& s, Eigen::Ref<Eigen::MatrixCmp> h, Eigen::Ref<Eigen::VectorCmp> g) const CASEMA_NOEXCEPT
{
	g.setZero();
	h.setZero();

	// Dini's expansion
	for (int n = 0; n < _besselZerosPhi.size(); ++n)
	{
		const mpfr::mpcomplex phiVal = phi(s) + _besselZerosPhi[n];
		const mpfr::mpcomplex alpha = sqrt(_one + _fourDaxOverUsq * phiVal);
		const mpfr::mpcomplex onePlusAlpha = _one + alpha;
		const mpfr::mpcomplex oneMinusAlpha = _one - alpha;
		const mpfr::mpcomplex outlet = _fourExpPeclet * alpha / (sqr(onePlusAlpha) * exp(_halfPeclet * onePlusAlpha) - sqr(oneMinusAlpha) * exp(_halfPeclet * oneMinusAlpha));

		for (int outZone = 0; outZone < _radialEdges.size()-1; ++outZone)
		{
			for (int inZone = 0; inZone < _radialEdges.size()-1; ++inZone)
				h(outZone, inZone) += outlet * _besselCoeffs(outZone, n) * _besselCoeffs(inZone, n) * _j0nSq[n];
		}
	}

	for (int outZone = 0; outZone < _radialEdges.size()-1; ++outZone)
		h.row(outZone) *= _crossAreaZone[outZone];
}

void GeneralRateModel2D::evaluate(const mpfr::mpcomplex& s, const mpfr::mpreal& z, Eigen::Ref<Eigen::MatrixCmp> h, Eigen::Ref<Eigen::VectorCmp> g) const CASEMA_NOEXCEPT
{
	evaluate(s, h, g);
}


void registerGeneralRateModel2D(std::unordered_map<std::string, std::function<UnitOperation*(int)>>& models)
{
	models[GeneralRateModel2D::identifier()] = [](int uoId) { return new GeneralRateModel2D(uoId); };
	models["GRM2D"] = [](int uoId) { return new GeneralRateModel2D(uoId); };
}

}  // namespace model

}  // namespace casema
