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

#include "model/ColumnLikeModel.hpp"
#include "io/ParameterProvider.hpp"
#include "model/ParamReaderScopes.hpp"
#include "model/ParamReaderUtils.hpp"
#include "Exceptions.hpp"

#include <numeric>

namespace casema
{

namespace model
{

/**
 * @brief Creates a ConvectionDispersionLikeModel
 */
ColumnLikeModel::ColumnLikeModel(int unitOpIdx) : ConvectionDispersionLikeModel(unitOpIdx)
{
}

ColumnLikeModel::~ColumnLikeModel() CASEMA_NOEXCEPT
{
}

/**
 * @brief Reads model parameters
 * @details Only reads parameters that do not affect model structure (e.g., discretization).
 * @param [in] paramProvider Parameter provider for reading parameters
 * @return @c true if configuration went fine, @c false otherwise
 */
bool ColumnLikeModel::configure(io::IParameterProvider& paramProvider)
{
	ConvectionDispersionLikeModel::configure(paramProvider);

	// Read geometry parameters
	const std::string unitType = paramProvider.getString("UNIT_TYPE");
	if (unitType != "DPFR")
	{
		if ((unitType != "LUMPED_RATE_MODEL_WITHOUT_PORES") && (unitType != "LRM"))
		{
			if ((unitType == "GENERAL_RATE_MODEL_2D") && paramProvider.isArray("COL_POROSITY"))
			{
				const std::vector<double> cp = paramProvider.getDoubleArray("COL_POROSITY");
				for (int i = 1; i < cp.size(); ++i)
				{
					if (cp[0] != cp[i])
						throw InvalidParameterException("Field COL_POROSITY has to have the same value if it is an array");
				}

				_colPorosity = cp[0];
			}
			else
			{
				if (paramProvider.isArray("COL_POROSITY"))
					throw InvalidParameterException("Only scalar COL_POROSITY supported");

				_colPorosity = paramProvider.getDouble("COL_POROSITY");
			}
		}
		else
			_colPorosity = paramProvider.getDouble("TOTAL_POROSITY");
	}
	else
		_colPorosity = 1.0;

	int numParType = 1;
	if (paramProvider.exists("PAR_TYPE_VOLFRAC"))
		numParType = paramProvider.numElements("PAR_TYPE_VOLFRAC");

	const std::vector<std::string> bindModelNames = paramProvider.getStringArray("ADSORPTION_MODEL");

	if ((bindModelNames.size() != 1) && (bindModelNames.size() != numParType))
		throw InvalidParameterException("Field ADSORPTION_MODEL contains too few elements (" + std::to_string(numParType) + " required)");

	_kA.reserve(numParType);
	_kD.reserve(numParType);
	_isKinetic.reserve(numParType);

	for (int i = 0; i < numParType; ++i)
	{
		if (bindModelNames[i] == "LINEAR")
		{
			MultiplexedScopeSelector scopeGuard(paramProvider, "adsorption", numParType == 1, i, true, true);

			if (paramProvider.numElements("LIN_KA") != 1)
				throw InvalidParameterException("Field LIN_KA must be scalar");
			if (paramProvider.numElements("LIN_KD") != 1)
				throw InvalidParameterException("Field LIN_KD must be scalar");
			if (paramProvider.numElements("IS_KINETIC") != 1)
				throw InvalidParameterException("Field IS_KINETIC must be scalar");

			_kA.push_back(paramProvider.getDouble("LIN_KA"));
			_kD.push_back(paramProvider.getDouble("LIN_KD"));
			_isKinetic.push_back(paramProvider.getBool("IS_KINETIC"));
		}
		else if (bindModelNames[i] == "NONE")
		{
			_kA.push_back(0.0);
			_kD.push_back(1.0);
			_isKinetic.push_back(false);
		}
		else
			throw InvalidParameterException("Field ADSORPTION_MODEL contains unsupported binding model (" + bindModelNames[i] + ")");
	}

	return true;
}

void ColumnLikeModel::setVelocityFromFlowRate(mpfr::mpreal const* in)
{
	if (_crossSection > 0.0)
		_velocity = in[0] / (_crossSection * _colPorosity);
}


/**
 * @brief Creates a ConvectionDispersionLikeModel
 */
ColumnWithParticles::ColumnWithParticles(int unitOpIdx) : ColumnLikeModel(unitOpIdx)
{
}

ColumnWithParticles::~ColumnWithParticles() CASEMA_NOEXCEPT
{
}

/**
 * @brief Reads model parameters
 * @details Only reads parameters that do not affect model structure (e.g., discretization).
 * @param [in] paramProvider Parameter provider for reading parameters
 * @return @c true if configuration went fine, @c false otherwise
 */
bool ColumnWithParticles::configure(io::IParameterProvider& paramProvider)
{
	ColumnLikeModel::configure(paramProvider);

	int numParType = 1;
	if (paramProvider.exists("PAR_TYPE_VOLFRAC"))
	{
		const std::vector<double> parTypeVolFrac = paramProvider.getDoubleArray("PAR_TYPE_VOLFRAC");
		numParType = paramProvider.numElements("PAR_TYPE_VOLFRAC");
		_parTypeVolFrac = std::move(toMPreal(parTypeVolFrac, numParType));
	}
	else
	{
		_parTypeVolFrac.push_back(1.0);
	}

	const std::vector<double> parRadius = paramProvider.getDoubleArray("PAR_RADIUS");
	if ((parRadius.size() != 1) && (parRadius.size() != numParType))
		throw InvalidParameterException("Field PAR_RADIUS has to have 1 or " + std::to_string(numParType) + " elements");
	_parRadius = std::move(toMPreal(parRadius, numParType));

	const std::vector<double> parPorosity = paramProvider.getDoubleArray("PAR_POROSITY");
	if ((parPorosity.size() != 1) && (parPorosity.size() != numParType))
		throw InvalidParameterException("Field PAR_POROSITY has to have 1 or " + std::to_string(numParType) + " elements");
	_parPorosity = std::move(toMPreal(parPorosity, numParType));

	const std::vector<double> filmDiffusion = paramProvider.getDoubleArray("FILM_DIFFUSION");
	if ((filmDiffusion.size() != 1) && (filmDiffusion.size() != numParType))
		throw InvalidParameterException("Field FILM_DIFFUSION has to have 1 or " + std::to_string(numParType) + " elements");
	_filmDiffusion = std::move(toMPreal(filmDiffusion, numParType));

	_betaC = (_one - _colPorosity) / _colPorosity;

	return true;
}

template <typename eval_t>
eval_t ColumnWithParticles::phiImpl(const eval_t& s) const CASEMA_NOEXCEPT
{
	eval_t p(0);
	for (int j = 0; j < _parTypeVolFrac.size(); ++j)
		p += _parTypeVolFrac[j] / _parRadius[j] * _filmDiffusion[j] * (_one - particleF(j, s));

	p *= 3 * _betaC;

	return p + s;
}

/**
 * @brief Creates a ConvectionDispersionLikeModel
 */
ColumnWithPoreDiffusion::ColumnWithPoreDiffusion(int unitOpIdx) : ColumnWithParticles(unitOpIdx)
{
}

ColumnWithPoreDiffusion::~ColumnWithPoreDiffusion() CASEMA_NOEXCEPT
{
}

/**
 * @brief Reads model parameters
 * @details Only reads parameters that do not affect model structure (e.g., discretization).
 * @param [in] paramProvider Parameter provider for reading parameters
 * @return @c true if configuration went fine, @c false otherwise
 */
bool ColumnWithPoreDiffusion::configure(io::IParameterProvider& paramProvider)
{
	ColumnWithParticles::configure(paramProvider);

	paramProvider.pushScope("discretization");
	_nBound = std::move(paramProvider.getIntArray("NBOUND"));
	paramProvider.popScope();

	const int numBound = std::accumulate(_nBound.begin(), _nBound.end(), 0);
	const int numParType = _parTypeVolFrac.size();

	if (paramProvider.exists("PAR_CORERADIUS"))
	{
		const std::vector<double> parCoreRadius = paramProvider.getDoubleArray("PAR_CORERADIUS");
		if ((parCoreRadius.size() != 1) && (parCoreRadius.size() != numParType))
			throw InvalidParameterException("Field PAR_CORERADIUS has to have 1 or " + std::to_string(numParType) + " elements");
		_parCoreRadius = std::move(toMPreal(parCoreRadius, numParType));
	}
	else
		_parCoreRadius = std::move(std::vector<mpfr::mpreal>(numParType, 0.0));

	const std::vector<double> parDiffusion = paramProvider.getDoubleArray("PAR_DIFFUSION");
	if ((parDiffusion.size() != 1) && (parDiffusion.size() != numParType))
		throw InvalidParameterException("Field PAR_DIFFUSION has to have 1 or " + std::to_string(numParType) + " elements");
	_parDiffusion = std::move(toMPreal(parDiffusion, numParType));

	const std::vector<double> parSurfDiffusion = paramProvider.getDoubleArray("PAR_SURFDIFFUSION");
	if ((parSurfDiffusion.size() != 1) && (parSurfDiffusion.size() != numBound))
		throw InvalidParameterException("Field PAR_SURFDIFFUSION has to have 1 or " + std::to_string(numBound) + " elements");
	_parSurfDiffusion = std::move(toMPreal(parSurfDiffusion, numBound));

	if (numBound != numParType)
	{
		// Introduce zero padding
		std::vector<mpfr::mpreal> p(numParType);
		int j = 0;
		for (int i = 0; i < numParType; ++i)
		{
			if (_nBound[i] > 0)
			{
				p[i] = _parSurfDiffusion[j];
				++j;
			}
			else
				p[i] = 0;
		}
		_parSurfDiffusion = std::move(p);
	}

	_temp1.reserve(_parTypeVolFrac.size());
	_temp2.reserve(_parTypeVolFrac.size());
	for (int j = 0; j < _parTypeVolFrac.size(); ++j)
	{
		if (_isKinetic[j])
		{
			if ((_parSurfDiffusion[j] == 0) || (_kA[j] == 0))
				_temp1.push_back(_parPorosity[j] * _parDiffusion[j] / (_filmDiffusion[j] * _parRadius[j]));
			else
				_temp1.push_back(mpfr::mpreal(4) * _parDiffusion[j] * _parSurfDiffusion[j]);

			_temp2.push_back(_kA[j] * (_one - _parPorosity[j]) / _parPorosity[j]);
		}
		else
		{
			const mpfr::mpreal W = _parPorosity[j] * _parDiffusion[j] + (_one - _parPorosity[j]) * _parSurfDiffusion[j] * _kA[j] / _kD[j];
			_temp1.push_back(W / (_filmDiffusion[j] * _parRadius[j]));
			_temp2.push_back(W / (_parPorosity[j] + (_one - _parPorosity[j]) * _kA[j] / _kD[j]));
		}
	}

	return true;
}

template <typename eval_t>
eval_t ColumnWithPoreDiffusion::particleFimpl(int j, const eval_t& s) const CASEMA_NOEXCEPT
{
	if (_isKinetic[j])
	{
		if ((_parSurfDiffusion[j] == 0) || (_kA[j] == 0))
		{
			const eval_t rpOverPsi = _parRadius[j] / sqrt((_parDiffusion[j] / s) * (s + _kD[j]) / (s + _kD[j] + _temp2[j]));
			return 1 / (_one + _temp1[j] * (rpOverPsi * coth(rpOverPsi) - _one));
		}
		else
		{
			const eval_t kappaLeft = _parSurfDiffusion[j] * (s + _temp2[j]);
			const eval_t kappaRight = _parDiffusion[j] * (s + _kD[j]);
			const eval_t kappa1 = kappaLeft + kappaRight;
			const eval_t kappa2 = kappaLeft - kappaRight;
			const eval_t gamma = sqr(kappa1) - _temp1[j] * s * (_temp2[j] + _kD[j] + s);
			const eval_t sqrtGamma = sqrt(gamma);
			const mpfr::mpreal denomV = 2 * _parDiffusion[j] * _kA[j];
			const eval_t v11 = -(kappa2 + sqrtGamma) / denomV;
			const eval_t v21 = -(kappa2 - sqrtGamma) / denomV;
			const eval_t denomTau = 2 * s * (_temp2[j] + _kD[j] + s);
			const eval_t tau1 = (kappa1 - sqrtGamma) / denomTau;
			const eval_t tau2 = (kappa1 + sqrtGamma) / denomTau;

			const eval_t sqrtTau1 = sqrt(tau1);
			const eval_t sqrtTau2 = sqrt(tau2);
			const eval_t rpOverSqrtTau1 = _parRadius[j] / sqrtTau1;
			const eval_t rpOverSqrtTau2 = _parRadius[j] / sqrtTau2;
			const eval_t coshRpOverSqrtTau1 = cosh(rpOverSqrtTau1);
			const eval_t coshRpOverSqrtTau2 = cosh(rpOverSqrtTau2);
			const eval_t sinhRpOverSqrtTau1 = sinh(rpOverSqrtTau1);
			const eval_t sinhRpOverSqrtTau2 = sinh(rpOverSqrtTau2);
			const eval_t nu1 = _parRadius[j] * coshRpOverSqrtTau1 - sqrtTau1 * sinhRpOverSqrtTau1;
			const eval_t nu2 = _parRadius[j] * coshRpOverSqrtTau2 - sqrtTau2 * sinhRpOverSqrtTau2;
			const eval_t nuTilde = sqrtTau1 * coshRpOverSqrtTau2 * sinhRpOverSqrtTau1 * v11 - sqrtTau2 * coshRpOverSqrtTau1 * sinhRpOverSqrtTau2 * v21;
			const eval_t mu = _parRadius[j] * nuTilde + sqrtGamma / _kA[j] * (sqrtTau1 * sqrtTau2 / _parDiffusion[j] * sinhRpOverSqrtTau1 * sinhRpOverSqrtTau2 - _parPorosity[j] / (_filmDiffusion[j] * _parRadius[j]) * nu1 * nu2);

			return (v11 * sqrtTau1 * sinhRpOverSqrtTau1 * nu2 - v21 * sqrtTau2 * sinhRpOverSqrtTau2 * nu1) / mu;
		}
	}
	else
	{
		const eval_t rpOverPsi = _parRadius[j] / sqrt(_temp2[j] / s);
		return 1 / (_one + _temp1[j] * (rpOverPsi * coth(rpOverPsi) - _one));
	}
}

}  // namespace model

}  // namespace casema
