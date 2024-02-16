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

/**
 * @file 
 * Defines the convection dispersion transport operator.
 */

#ifndef CASEMA_COLUMNLIKE_HPP_
#define CASEMA_COLUMNLIKE_HPP_

#include "model/ConvectionDispersionLikeModel.hpp"

namespace casema
{

namespace model
{

class ColumnLikeModel : public ConvectionDispersionLikeModel
{
public:

	ColumnLikeModel(int unitOpIdx);
	~ColumnLikeModel() CASEMA_NOEXCEPT;

	virtual bool configure(io::IParameterProvider& paramProvider);

protected:

	virtual void setVelocityFromFlowRate(mpfr::mpreal const* in);

	mpfr::mpreal _colPorosity; //!< Column porosity (external porosity) \f$ \varepsilon_c \f$
	std::vector<mpfr::mpreal> _kA; //!< Adsorption rate
	std::vector<mpfr::mpreal> _kD; //!< Desorption rate
	std::vector<bool> _isKinetic; //!< Desorption rate
};

class ColumnWithParticles : public ColumnLikeModel
{
public:

	ColumnWithParticles(int unitOpIdx);
	~ColumnWithParticles() CASEMA_NOEXCEPT;

	virtual bool configure(io::IParameterProvider& paramProvider);

protected:

	virtual mpfr::mpreal particleF(int j, const mpfr::mpreal& s) const CASEMA_NOEXCEPT = 0;
	virtual mpfr::mpcomplex particleF(int j, const mpfr::mpcomplex& s) const CASEMA_NOEXCEPT = 0;

	virtual mpfr::mpreal phi(const mpfr::mpreal& s) const CASEMA_NOEXCEPT { return phiImpl(s); }
	virtual mpfr::mpcomplex phi(const mpfr::mpcomplex& s) const CASEMA_NOEXCEPT { return phiImpl(s); }

	std::vector<mpfr::mpreal> _parRadius; //!< Particle radius \f$ r_p \f$
	std::vector<mpfr::mpreal> _parPorosity; //!< Particle porosity (internal porosity) \f$ \varepsilon_p \f$
	std::vector<mpfr::mpreal> _parTypeVolFrac; //!< Volume fraction of each particle type

	// Vectorial parameters
	std::vector<mpfr::mpreal> _filmDiffusion; //!< Film diffusion coefficient \f$ k_f \f$

	mpfr::mpreal _betaC;

private:

	template <typename eval_t>
	eval_t phiImpl(const eval_t& s) const CASEMA_NOEXCEPT;
};

class ColumnWithPoreDiffusion : public ColumnWithParticles
{
public:

	ColumnWithPoreDiffusion(int unitOpIdx);
	~ColumnWithPoreDiffusion() CASEMA_NOEXCEPT;

	virtual bool configure(io::IParameterProvider& paramProvider);

protected:

	virtual mpfr::mpreal particleF(int j, const mpfr::mpreal& s) const CASEMA_NOEXCEPT { return particleFimpl(j, s); }
	virtual mpfr::mpcomplex particleF(int j, const mpfr::mpcomplex& s) const CASEMA_NOEXCEPT { return particleFimpl(j, s); }

	std::vector<mpfr::mpreal> _parCoreRadius; //!< Particle core radius \f$ r_c \f$
	std::vector<mpfr::mpreal> _parDiffusion; //!< Particle diffusion coefficient \f$ D_p \f$
	std::vector<mpfr::mpreal> _parSurfDiffusion; //!< Particle surface diffusion coefficient \f$ D_s \f$
	std::vector<int> _nBound;

	std::vector<mpfr::mpreal> _temp1; //!< 
	std::vector<mpfr::mpreal> _temp2; //!< 

private:

	template <typename eval_t>
	eval_t particleFimpl(int j, const eval_t& s) const CASEMA_NOEXCEPT;
};

} // namespace model
} // namespace casema

#endif  // CASEMA_COLUMNLIKE_HPP_
