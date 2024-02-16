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
 * Defines an unit operation model interface.
 */

#ifndef CASEMA_UNITOPERATION_HPP_
#define CASEMA_UNITOPERATION_HPP_

#include "casemaCompilerInfo.hpp"
#include "MPReal.hpp"
#include "MPComplex.hpp"
#include "MPComplexEigen.hpp"

#include <tuple>
#include <vector>

namespace casema
{

namespace io
{
	class IParameterProvider;
}

namespace model
{

/**
 * @brief Defines an unit operation model interface
 * @details 
 */
class UnitOperation
{
public:

	UnitOperation(int idxUnit) : _unitOpIdx(idxUnit) { }
	virtual ~UnitOperation() CASEMA_NOEXCEPT { }

	int unitOperationId() const CASEMA_NOEXCEPT { return _unitOpIdx; }	

	/**
	 * @brief (Re-)configures the model by extracting all non-structural parameters (e.g., model parameters) from the given @p paramProvider
	 * @details The scope of the casema::IParameterProvider is left unchanged on return.
	 *          
	 *          The structure of the model is left unchanged, that is, the number of degrees of
	 *          freedom stays the same. Only true (non-structural) model parameters are read and
	 *          changed. Parameters that concern discretization (e.g., number of cells), model
	 *          structure (e.g., number of components, binding model), and numerical solution
	 *          (e.g., tolerances in GMRES iterations) are left untouched.
	 *          
	 *          This function may only be called if configureModelDiscretization() has been called
	 *          in the past. Contrary to configureModelDiscretization(), it can be called multiple
	 *          times.
	 * 
	 * @param [in] paramProvider Parameter provider
	 * @return @c true if the configuration was successful, otherwise @c false
	 */
	virtual bool configure(io::IParameterProvider& paramProvider);

	/**
	 * @brief Returns the number of components
	 * @details It is assumed that the number of components is also the number of inputs
	 *          and outputs of the unit operation.
	 * @return Number of components
	 */
	virtual int numComponents() const CASEMA_NOEXCEPT { return _nComp; }

	virtual const char* unitOperationName() const CASEMA_NOEXCEPT = 0;

	/**
	 * @brief Returns whether this unit operation possesses an inlet
	 * @return @c true if the unit operation can take in a stream, otherwise @c false
	 */
	virtual bool hasInlet() const CASEMA_NOEXCEPT = 0;

	/**
	 * @brief Returns whether this unit operation possesses an outlet
	 * @return @c true if the unit operation can output a stream, otherwise @c false
	 */
	virtual bool hasOutlet() const CASEMA_NOEXCEPT = 0;

	/**
	 * @brief Returns the number of inlet ports
	 * @return Number of inlet ports
	 */
	virtual int numInletPorts() const CASEMA_NOEXCEPT = 0;

	/**
	 * @brief Returns the number of outlet ports
	 * @return Number of outlet ports
	 */
	virtual int numOutletPorts() const CASEMA_NOEXCEPT = 0;

	/**
	 * @brief Sets the flow rates for the current time section
	 * @details The flow rates may change due to valve switches.
	 * @param [in] in Array with total volumetric inlet flow rate for each port
	 * @param [in] out Array with total volumetric outlet flow rate for each port
	 */
	virtual void setFlowRates(mpfr::mpreal const* in, mpfr::mpreal const* out) = 0;

	virtual void setBesselZeros(mpfr::mpreal const* zeros, int n) CASEMA_NOEXCEPT { }
	virtual bool needsBesselZeros() const CASEMA_NOEXCEPT { return false; }

	/**
	* @brief Returns whether this unit operation supports non-matching volumetric inlet and outlet flow rates
	* @details If inlet and outlet flow rates do not match, mass is accumulated or lost during time integration.
	* @return @c true if flow rates are allowed to differ, otherwise @c false
	*/
	virtual bool canAccumulate() const CASEMA_NOEXCEPT = 0;

	virtual void setSectionTimes(mpfr::mpreal const* secTimes, int nTimes) CASEMA_NOEXCEPT = 0;

	virtual void evaluate(const mpfr::mpcomplex& s, Eigen::Ref<Eigen::MatrixCmp> h, Eigen::Ref<Eigen::VectorCmp> g) const CASEMA_NOEXCEPT = 0;
	virtual void evaluate(const mpfr::mpcomplex& s, const mpfr::mpreal& z, Eigen::Ref<Eigen::MatrixCmp> h, Eigen::Ref<Eigen::VectorCmp> g) const CASEMA_NOEXCEPT = 0;

	virtual bool hasValidEstimate() const CASEMA_NOEXCEPT = 0;
	virtual mpfr::mpreal estimate(const mpfr::mpreal& abscissa) const CASEMA_NOEXCEPT = 0;
	virtual mpfr::mpreal timeDomainUpperBound(mpfr::mpreal const* in) const CASEMA_NOEXCEPT = 0;
	virtual mpfr::mpreal truncationError(const mpfr::mpreal& M, const mpfr::mpreal& abscissa, const mpfr::mpreal& T, const mpfr::mpreal& nTerms) const CASEMA_NOEXCEPT = 0;
	virtual mpfr::mpreal inverseTruncationError(const mpfr::mpreal& M, const mpfr::mpreal& abscissa, const mpfr::mpreal& T, const mpfr::mpreal& error) const CASEMA_NOEXCEPT = 0;

protected:
	int _unitOpIdx; //!< Unit operation index
	int _nComp;
};

} // namespace model
} // namespace casema

#endif  // CASEMA_UNITOPERATION_HPP_
