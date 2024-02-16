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

#include "model/LumpedRateModelWithPores.hpp"
#include "io/ParameterProvider.hpp"
#include "Exceptions.hpp"

#include <unordered_map>

namespace casema
{

namespace model
{

/**
 * @brief Creates a ConvectionDispersionLikeModel
 */
LumpedRateModelWithPores::LumpedRateModelWithPores(int unitOpIdx) : ColumnWithParticles(unitOpIdx)
{
}

LumpedRateModelWithPores::~LumpedRateModelWithPores() CASEMA_NOEXCEPT
{
}

/**
 * @brief Reads model parameters
 * @details Only reads parameters that do not affect model structure (e.g., discretization).
 * @param [in] paramProvider Parameter provider for reading parameters
 * @return @c true if configuration went fine, @c false otherwise
 */
bool LumpedRateModelWithPores::configure(io::IParameterProvider& paramProvider)
{
	ColumnWithParticles::configure(paramProvider);

	_threeKfOverRp.reserve(_parTypeVolFrac.size());
	_parFactor.reserve(_parTypeVolFrac.size());
	for (int j = 0; j < _parTypeVolFrac.size(); ++j)
	{
		_threeKfOverRp.push_back(3 * _filmDiffusion[j] / _parRadius[j]);
		if (_isKinetic[j])
			_parFactor.push_back((_one - _parPorosity[j]) * _kA[j]);
		else
			_parFactor.push_back(_parPorosity[j] + (_one - _parPorosity[j]) * _kA[j] / _kD[j]);
	}

	return true;
}

template <typename eval_t>
eval_t LumpedRateModelWithPores::particleFimpl(int j, const eval_t& s) const CASEMA_NOEXCEPT
{
	if (_isKinetic[j])
		return _threeKfOverRp[j] / (_threeKfOverRp[j] + s * (_parPorosity[j] + _parFactor[j] / (_kD[j] + s)));
	else
		return _threeKfOverRp[j] / (_threeKfOverRp[j] + s * _parFactor[j]);
}

void registerLumpedRateModelWithPores(std::unordered_map<std::string, std::function<UnitOperation*(int)>>& models)
{
	models[LumpedRateModelWithPores::identifier()] = [](int uoId) { return new LumpedRateModelWithPores(uoId); };
	models["LRMP"] = [](int uoId) { return new LumpedRateModelWithPores(uoId); };
}

}  // namespace model

}  // namespace casema
