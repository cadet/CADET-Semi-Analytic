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

#include "model/StirredTankModel.hpp"
#include "io/ParameterProvider.hpp"
#include "Exceptions.hpp"
#include "Constants.hpp"

#include <unordered_map>
#include <functional>
#include <limits>
#include <string>

namespace casema
{

namespace model
{

CSTRModel::CSTRModel(int unitOpIdx) : UnitOperation(unitOpIdx)
{
}

CSTRModel::~CSTRModel() CASEMA_NOEXCEPT
{
}

void CSTRModel::setFlowRates(mpfr::mpreal const* in, mpfr::mpreal const* out) 
{ 
	_flowRateIn = in[0];
	_flowRateOut = out[0];

//	if (abs(_flowRateIn - (_flowRateOut + _flowRateFilter)) > std::numeric_limits<mpfr::mpreal>::epsilon())
//		throw InvalidParameterException("CSTR does not support variable volume");
}

bool CSTRModel::configure(io::IParameterProvider& paramProvider)
{
	UnitOperation::configure(paramProvider);

	_flowRateFilter = 0.0;
	if (paramProvider.exists("FLOWRATE_FILTER"))
	{
		if (paramProvider.isArray("FLOWRATE_FILTER"))
			throw InvalidParameterException("Field FLOWRATE_FILTER has to be scalar");
		_flowRateFilter = paramProvider.getDouble("FLOWRATE_FILTER");
	}

	_initC = paramProvider.getDouble("INIT_C");
	_initV = paramProvider.getDouble("INIT_VOLUME");

	if (_initV <= 0.0)
		throw InvalidParameterException("CSTR does not support zero volume");

	return true;
}

void CSTRModel::evaluate(const mpfr::mpcomplex& s, Eigen::Ref<Eigen::MatrixCmp> h, Eigen::Ref<Eigen::VectorCmp> g) const CASEMA_NOEXCEPT
{
	const mpfr::mpcomplex denom = _initV * s + _flowRateOut;
	h(0,0) = _flowRateIn / denom;
	g(0) = _initV * _initC / denom;
}

void CSTRModel::evaluate(const mpfr::mpcomplex& s, const mpfr::mpreal& z, Eigen::Ref<Eigen::MatrixCmp> h, Eigen::Ref<Eigen::VectorCmp> g) const CASEMA_NOEXCEPT
{
	evaluate(s, h, g);
}

bool CSTRModel::hasValidEstimate() const CASEMA_NOEXCEPT
{
	return (_initC == 0.0) && (_initV > 0.0) && (abs(_flowRateIn - (_flowRateOut + _flowRateFilter)) <= std::numeric_limits<mpfr::mpreal>::epsilon());
}

mpfr::mpreal CSTRModel::estimate(const mpfr::mpreal& abscissa) const CASEMA_NOEXCEPT
{
	return _flowRateIn / _initV;
}

mpfr::mpreal CSTRModel::truncationError(const mpfr::mpreal& M, const mpfr::mpreal& abscissa, const mpfr::mpreal& T, const mpfr::mpreal& nTerms) const CASEMA_NOEXCEPT
{
	const mpfr::mpreal qOutAV0 = _flowRateOut + _initV * abscissa;
	return 2 * M * T * _flowRateIn / (Constants<mpfr::mpreal>::pi() * qOutAV0) * asinh(qOutAV0 * T / (Constants<mpfr::mpreal>::pi() * _initV * nTerms));
}

mpfr::mpreal CSTRModel::inverseTruncationError(const mpfr::mpreal& M, const mpfr::mpreal& abscissa, const mpfr::mpreal& T, const mpfr::mpreal& error) const CASEMA_NOEXCEPT
{
	const mpfr::mpreal qOutAV0 = _flowRateOut + _initV * abscissa;
	return ceil( T * qOutAV0 / (_initV * Constants<mpfr::mpreal>::pi() * sinh( error * Constants<mpfr::mpreal>::pi() * qOutAV0 / (2 * M * T * _flowRateIn) )) );
}


void registerCSTRModel(std::unordered_map<std::string, std::function<UnitOperation*(int)>>& models)
{
	models[CSTRModel::identifier()] = [](int uoId) { return new CSTRModel(uoId); };
}

}  // namespace model

}  // namespace casema
