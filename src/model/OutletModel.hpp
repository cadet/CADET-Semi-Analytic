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
 * Defines the outlet model.
 */

#ifndef CASEMA_OUTLETMODEL_HPP_
#define CASEMA_OUTLETMODEL_HPP_

#include "model/UnitOperation.hpp"

namespace casema
{

namespace model
{

/**
 * @brief Outlet model
 * @details This unit operation model is solely for buffering concentration profiles
 *          for readout.
 */
class OutletModel : public UnitOperation
{
public:

	OutletModel(int unitOpIdx) : UnitOperation(unitOpIdx) { }
	virtual ~OutletModel() CASEMA_NOEXCEPT { }

	virtual int numInletPorts() const CASEMA_NOEXCEPT { return 1; }
	virtual int numOutletPorts() const CASEMA_NOEXCEPT { return 0; }
	virtual bool hasInlet() const CASEMA_NOEXCEPT { return true; }
	virtual bool hasOutlet() const CASEMA_NOEXCEPT { return false; }
	virtual bool canAccumulate() const CASEMA_NOEXCEPT { return true; }

	static const char* identifier() { return "OUTLET"; }
	virtual const char* unitOperationName() const CASEMA_NOEXCEPT { return "OUTLET"; }

	virtual bool configure(io::IParameterProvider& paramProvider);
	virtual void setFlowRates(mpfr::mpreal const* in, mpfr::mpreal const* out) { }
	virtual void setSectionTimes(mpfr::mpreal const* secTimes, int nTimes) CASEMA_NOEXCEPT { }

	virtual void evaluate(const mpfr::mpcomplex& s, Eigen::Ref<Eigen::MatrixCmp> h, Eigen::Ref<Eigen::VectorCmp> g) const CASEMA_NOEXCEPT { h(0,0) = 1.0; g(0) = 0.0; }
	virtual void evaluate(const mpfr::mpcomplex& s, const mpfr::mpreal& z, Eigen::Ref<Eigen::MatrixCmp> h, Eigen::Ref<Eigen::VectorCmp> g) const CASEMA_NOEXCEPT { h(0,0) = 1.0; g(0) = 0.0; }

	virtual bool hasValidEstimate() const CASEMA_NOEXCEPT { return true; }
	virtual mpfr::mpreal estimate(const mpfr::mpreal& abscissa) const CASEMA_NOEXCEPT { return 1.0; }
	virtual mpfr::mpreal timeDomainUpperBound(mpfr::mpreal const* in) const CASEMA_NOEXCEPT { return in[0]; }
	virtual mpfr::mpreal truncationError(const mpfr::mpreal& M, const mpfr::mpreal& abscissa, const mpfr::mpreal& T, const mpfr::mpreal& nTerms) const CASEMA_NOEXCEPT { return -1; }
	virtual mpfr::mpreal inverseTruncationError(const mpfr::mpreal& M, const mpfr::mpreal& abscissa, const mpfr::mpreal& T, const mpfr::mpreal& error) const CASEMA_NOEXCEPT { return -1; }
};

} // namespace model
} // namespace casema

#endif  // CASEMA_OUTLETMODEL_HPP_
