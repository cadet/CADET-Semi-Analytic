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

#include "model/ConvectionDispersionLikeModel.hpp"
#include "io/ParameterProvider.hpp"
#include "Exceptions.hpp"
#include "Constants.hpp"
#include "ExpInt.hpp"

namespace casema
{

namespace model
{

/**
 * @brief Creates a ConvectionDispersionLikeModel
 */
ConvectionDispersionLikeModel::ConvectionDispersionLikeModel(int unitOpIdx) : UnitOperation(unitOpIdx), _one("1.0")
{
}

ConvectionDispersionLikeModel::~ConvectionDispersionLikeModel() CASEMA_NOEXCEPT
{
}

/**
 * @brief Reads model parameters
 * @details Only reads parameters that do not affect model structure (e.g., discretization).
 * @param [in] int Unit operation id of the owning unit operation model
 * @param [in] paramProvider Parameter provider for reading parameters
 * @param [out] parameters Map in which local parameters are inserted
 * @return @c true if configuration went fine, @c false otherwise
 */
bool ConvectionDispersionLikeModel::configure(io::IParameterProvider& paramProvider)
{
	UnitOperation::configure(paramProvider);

	// Read geometry parameters
	_colLength = paramProvider.getDouble("COL_LENGTH");

	// Read cross section area or set to -1
	_crossSection = -1.0;
	if (paramProvider.exists("CROSS_SECTION_AREA"))
		_crossSection = paramProvider.getDouble("CROSS_SECTION_AREA");

	// Read VELOCITY
	_velocity = -1.0;
	if (paramProvider.exists("VELOCITY"))
	{
		if ((paramProvider.getString("UNIT_TYPE") == "GENERAL_RATE_MODEL_2D") && paramProvider.isArray("VELOCITY"))
		{
			const std::vector<double> v = paramProvider.getDoubleArray("VELOCITY");
			for (int i = 1; i < v.size(); ++i)
			{
				if (v[0] != v[i])
					throw InvalidParameterException("Field VELOCITY has to have the same value if it is an array");
			}

			_velocity = v[0];
		}
		else
		{
			if (paramProvider.isArray("VELOCITY"))
				throw InvalidParameterException("Only scalar VELOCITY supported");

			_velocity = paramProvider.getDouble("VELOCITY");
		}
	}

	_colDispersion = -1.0;
	if (paramProvider.exists("COL_DISPERSION"))
	{
		if (paramProvider.isArray("COL_DISPERSION"))
			throw InvalidParameterException("Only scalar COL_DISPERSION supported");
		
		_colDispersion = paramProvider.getDouble("COL_DISPERSION");
	}
	else
		throw InvalidParameterException("Parameter COL_DISPERSION is required");

/*
	if ((_velocity < 0) && (_crossSection <= 0.0))
	{
		throw InvalidParameterException("At least one of CROSS_SECTION_AREA and VELOCITY has to be set");
	}
*/

	_initC = paramProvider.getDouble("INIT_C");

	return true;
}

void ConvectionDispersionLikeModel::setFlowRates(mpfr::mpreal const* in, mpfr::mpreal const* out)
{
	setVelocityFromFlowRate(in);
	onUpdatedVelocity();
}

void ConvectionDispersionLikeModel::onUpdatedVelocity() CASEMA_NOEXCEPT
{
	_fourExpPeclet = 4 * exp(_colLength * _velocity / _colDispersion);
	_fourDaxOverUsq = 4 * _colDispersion / sqr(_velocity);
	_halfPeclet = _colLength * _velocity / (_colDispersion * 2);
}

void ConvectionDispersionLikeModel::setVelocityFromFlowRate(mpfr::mpreal const* in)
{
	if (_crossSection > 0.0)
		_velocity = in[0] / _crossSection;
}

void ConvectionDispersionLikeModel::evaluate(const mpfr::mpcomplex& s, Eigen::Ref<Eigen::MatrixCmp> h, Eigen::Ref<Eigen::VectorCmp> g) const CASEMA_NOEXCEPT
{
	const mpfr::mpcomplex phiVal = phi(s);
	const mpfr::mpcomplex alpha = sqrt(_one + _fourDaxOverUsq * phiVal);
	const mpfr::mpcomplex onePlusAlpha = _one + alpha;
	const mpfr::mpcomplex oneMinusAlpha = _one - alpha;
	
	h(0,0) = _fourExpPeclet * alpha / (sqr(onePlusAlpha) * exp(_halfPeclet * onePlusAlpha) - sqr(oneMinusAlpha) * exp(_halfPeclet * oneMinusAlpha));
	g(0) = (_one - h(0,0)) * _initC / phiVal;
}

void ConvectionDispersionLikeModel::evaluate(const mpfr::mpcomplex& s, const mpfr::mpreal& z, Eigen::Ref<Eigen::MatrixCmp> h, Eigen::Ref<Eigen::VectorCmp> g) const CASEMA_NOEXCEPT
{
	evaluate(s, h, g);
}

bool ConvectionDispersionLikeModel::hasValidEstimate() const CASEMA_NOEXCEPT
{
	return (_initC == 0.0);
}

mpfr::mpreal ConvectionDispersionLikeModel::estimate(const mpfr::mpreal& abscissa) const CASEMA_NOEXCEPT
{
	const mpfr::mpreal two(2);
	const mpfr::mpreal halfPeclet = _colLength * _velocity / (two * _colDispersion);
	return mpfr::mpreal(8) * sqrt(two) * exp(halfPeclet - two) * _colDispersion / sqr(_colLength);
}

mpfr::mpreal ConvectionDispersionLikeModel::truncationError(const mpfr::mpreal& M, const mpfr::mpreal& abscissa, const mpfr::mpreal& T, const mpfr::mpreal& nTerms) const CASEMA_NOEXCEPT
{
	const mpfr::mpreal gamma = _colLength * sqrt(Constants<mpfr::mpreal>::pi() / (2 * T * _colDispersion));
	return sqrt(mpfr::mpreal("32")) * M * T / Constants<mpfr::mpreal>::pi() * exp(_halfPeclet) * expInt(gamma * sqrt(nTerms));
}

mpfr::mpreal ConvectionDispersionLikeModel::inverseTruncationError(const mpfr::mpreal& M, const mpfr::mpreal& abscissa, const mpfr::mpreal& T, const mpfr::mpreal& error) const CASEMA_NOEXCEPT
{
	const mpfr::mpreal gamma = _colLength * sqrt(Constants<mpfr::mpreal>::pi() / (2 * T * _colDispersion));
	return ceil( sqr( inverseExpInt(error * Constants<mpfr::mpreal>::pi() / ( M * T * sqrt(mpfr::mpreal("32")) * exp(_halfPeclet) )) / gamma ) );
}

}  // namespace model

}  // namespace casema
