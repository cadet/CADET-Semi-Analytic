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

#ifndef CASEMA_CONVECTIONDISPERSIONLIKEMODEL_HPP_
#define CASEMA_CONVECTIONDISPERSIONLIKEMODEL_HPP_

#include "model/UnitOperation.hpp"

#include <vector>

namespace casema
{

namespace model
{

class ConvectionDispersionLikeModel : public UnitOperation
{
public:

	ConvectionDispersionLikeModel(int unitOpIdx);
	~ConvectionDispersionLikeModel() CASEMA_NOEXCEPT;

	virtual bool configure(io::IParameterProvider& paramProvider);
	virtual void setFlowRates(mpfr::mpreal const* in, mpfr::mpreal const* out);
	virtual void setSectionTimes(mpfr::mpreal const* secTimes, int nTimes) CASEMA_NOEXCEPT { }

	virtual void evaluate(const mpfr::mpcomplex& s, Eigen::Ref<Eigen::MatrixCmp> h, Eigen::Ref<Eigen::VectorCmp> g) const CASEMA_NOEXCEPT;
	virtual void evaluate(const mpfr::mpcomplex& s, const mpfr::mpreal& z, Eigen::Ref<Eigen::MatrixCmp> h, Eigen::Ref<Eigen::VectorCmp> g) const CASEMA_NOEXCEPT;

	virtual bool hasValidEstimate() const CASEMA_NOEXCEPT;
	virtual mpfr::mpreal estimate(const mpfr::mpreal& abscissa) const CASEMA_NOEXCEPT;
	virtual mpfr::mpreal timeDomainUpperBound(mpfr::mpreal const* in) const CASEMA_NOEXCEPT { return in[0] + _initC; }
	virtual mpfr::mpreal truncationError(const mpfr::mpreal& M, const mpfr::mpreal& abscissa, const mpfr::mpreal& T, const mpfr::mpreal& nTerms) const CASEMA_NOEXCEPT;
	virtual mpfr::mpreal inverseTruncationError(const mpfr::mpreal& M, const mpfr::mpreal& abscissa, const mpfr::mpreal& T, const mpfr::mpreal& error) const CASEMA_NOEXCEPT;

	void velocity(const mpfr::mpreal& velocity) CASEMA_NOEXCEPT { _velocity = velocity; onUpdatedVelocity(); }
	mpfr::mpreal velocity() const CASEMA_NOEXCEPT { return _velocity; }

protected:

	void onUpdatedVelocity() CASEMA_NOEXCEPT;

	virtual mpfr::mpreal phi(const mpfr::mpreal& s) const CASEMA_NOEXCEPT { return s; }
	virtual mpfr::mpcomplex phi(const mpfr::mpcomplex& s) const CASEMA_NOEXCEPT { return s; }
	virtual void setVelocityFromFlowRate(mpfr::mpreal const* in);

	mpfr::mpreal _colLength; //!< Column length \f$ L \f$
	mpfr::mpreal _crossSection; //!< Cross section area 

	// Section dependent parameters
	mpfr::mpreal _colDispersion; //!< Column dispersion (may be section dependent) \f$ D_{\text{ax}} \f$
	mpfr::mpreal _velocity; //!< Interstitial velocity (may be section dependent) \f$ u \f$
	mpfr::mpreal _curVelocity; //!< Current interstitial velocity \f$ u \f$ in this time section

	mpfr::mpreal _initC; //!< Liquid bulk phase initial conditions

	mpfr::mpreal _fourExpPeclet; //!< Liquid bulk phase initial conditions
	mpfr::mpreal _fourDaxOverUsq; //!< Liquid bulk phase initial conditions
	mpfr::mpreal _halfPeclet; //!< Liquid bulk phase initial conditions
	const mpfr::mpreal _one; //!< Liquid bulk phase initial conditions
};

} // namespace model
} // namespace casema

#endif  // CASEMA_CONVECTIONDISPERSIONLIKEMODEL_HPP_
