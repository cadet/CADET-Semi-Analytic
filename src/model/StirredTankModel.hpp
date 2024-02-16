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
* Defines the CSTR model as a unit operation.
*/

#ifndef CASEMA_CSTR_HPP_
#define CASEMA_CSTR_HPP_

#include "model/UnitOperation.hpp"

#include <vector>

namespace casema
{

namespace model
{

/**
 * @brief Continuous stirred tank (reactor) model
 * @details This is a simple CSTR model with variable volume using the ``well mixed assumption''.
 * @f[\begin{align}
	\frac{\mathrm{d}}{\mathrm{d} t}\left( V \left[ c_i + \frac{1}{\beta} \sum_{j=1}^{N_{\text{bnd},i}} q_{i,j} \right] \right) &= F_{\text{in}} c_{\text{in},i} - F_{\text{out}} c_i \\
	a \frac{\partial q_{i,j}}{\partial t} &= f_{\text{iso},i,j}(c, q) \\
	\frac{\partial V}{\partial t} &= F_{\text{in}} - F_{\text{out}} - F_{\text{filter}}
\end{align} @f]
 * The model can be used as a plain stir tank without any binding states.
 */
class CSTRModel : public UnitOperation
{
public:

	CSTRModel(int unitOpIdx);
	virtual ~CSTRModel() CASEMA_NOEXCEPT;

	virtual int numInletPorts() const CASEMA_NOEXCEPT { return 1; }
	virtual int numOutletPorts() const CASEMA_NOEXCEPT { return 1; }
	virtual bool hasInlet() const CASEMA_NOEXCEPT { return true; }
	virtual bool hasOutlet() const CASEMA_NOEXCEPT { return true; }
	virtual bool canAccumulate() const CASEMA_NOEXCEPT { return true; }

	static const char* identifier() { return "CSTR"; }
	virtual const char* unitOperationName() const CASEMA_NOEXCEPT { return "CSTR"; }

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

protected:

	mpfr::mpreal _flowRateIn; //!< Volumetric flow rate of incoming stream
	mpfr::mpreal _flowRateOut; //!< Volumetric flow rate of drawn outgoing stream
	mpfr::mpreal _flowRateFilter; //!< Volumetric flow rate of liquid outtake stream

	mpfr::mpreal _initC; //!< Initial conditions, ordering: Liquid phase concentration, solid phase concentration, volume
	mpfr::mpreal _initV; //!< Initial conditions, ordering: Liquid phase concentration, solid phase concentration, volume

};

} // namespace model
} // namespace casema

#endif  // CASEMA_CSTR_HPP_
