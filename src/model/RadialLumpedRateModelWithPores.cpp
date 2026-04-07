// =============================================================================
//  CADET-semi-analytic - The semi-analytic extension of CADET
//
//  Copyright © 2015-present: Samuel Leweke¹² and the CADET-Semi-Analytic
//  authors, see the AUTHORS file
//
//    ¹ Forschungszentrum Juelich GmbH, IBG-1, Juelich, Germany.
//    ² University of Cologne, Cologne, Germany.
//
//  All rights reserved. This program and the accompanying materials
//  are made available under the terms of the GNU Public License v3.0 (or, at
//  your option, any later version) which accompanies this distribution, and
//  is available at http://www.gnu.org/licenses/gpl.html
// =============================================================================

#include "model/RadialLumpedRateModelWithPores.hpp"
#include "io/ParameterProvider.hpp"
#include "Exceptions.hpp"
#include "Constants.hpp"

#include <unordered_map>

namespace casema
{

namespace model
{

/**
 * @brief Creates a RadialLumpedRateModelWithPores
 */
RadialLumpedRateModelWithPores::RadialLumpedRateModelWithPores(int unitOpIdx)
	: LumpedRateModelWithPores(unitOpIdx), _flowDirection(1)
{
}

RadialLumpedRateModelWithPores::~RadialLumpedRateModelWithPores() CASEMA_NOEXCEPT
{
}

/**
 * @brief Reads model parameters
 * @details Only reads parameters that do not affect model structure (e.g., discretization).
 * @param [in] paramProvider Parameter provider for reading parameters
 * @return @c true if configuration went fine, @c false otherwise
 */
bool RadialLumpedRateModelWithPores::configure(io::IParameterProvider& paramProvider)
{
	LumpedRateModelWithPores::configure(paramProvider);

	_innerRadius = paramProvider.getDouble("COL_RADIUS_INNER");
	_outerRadius = paramProvider.getDouble("COL_RADIUS_OUTER");

	if (_innerRadius <= 0.0 || _outerRadius <= 0.0)
		throw InvalidParameterException("COL_RADIUS_INNER and COL_RADIUS_OUTER must be > 0");
	if (_innerRadius >= _outerRadius)
		throw InvalidParameterException("COL_RADIUS_INNER must be < COL_RADIUS_OUTER");

	_physColLength = _colLength;

	_flowDirection = 1;
	if (paramProvider.exists("VELOCITY_COEFF"))
	{
		const double velCoeff = paramProvider.getDouble("VELOCITY_COEFF");
		_flowDirection = (velCoeff >= 0.0) ? 1 : -1;
	}

	_colLength = _outerRadius - _innerRadius;

	return true;
}

void RadialLumpedRateModelWithPores::setVelocityFromFlowRate(mpfr::mpreal const* in)
{
	const mpfr::mpreal twoPi = 2 * Constants<mpfr::mpreal>::pi();
	const mpfr::mpreal rMean = (_innerRadius + _outerRadius) / 2;
	_velocity = mpfr::mpreal(_flowDirection) * in[0] / (twoPi * _physColLength * _colPorosity * rMean);
}


void registerRadialLumpedRateModelWithPores(std::unordered_map<std::string, std::function<UnitOperation*(int)>>& models)
{
	const std::function<UnitOperation*(int)> f = [](int uoId) { return new RadialLumpedRateModelWithPores(uoId); };
	models[RadialLumpedRateModelWithPores::identifier()] = f;
	models["rLRMP"] = f;
}

}  // namespace model

}  // namespace casema