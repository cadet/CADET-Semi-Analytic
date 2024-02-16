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
 * Defines the inlet model which encapsulates an inlet profile as unit operation.
 */

#ifndef CASEMA_INLETMODEL_HPP_
#define CASEMA_INLETMODEL_HPP_

#include "model/UnitOperation.hpp"

#include <vector>

namespace casema
{

namespace model
{

/**
 * @brief Inlet model which encapsulates an IInletProfile as unit operation
 * @details This unit operation model does not possess DOFs. Its values (inlet profile)
 *          is directly injected into other models' DOFs.
 */
class InletModel : public UnitOperation
{
public:

	InletModel(int unitOpIdx);
	virtual ~InletModel() CASEMA_NOEXCEPT;

	virtual int numInletPorts() const CASEMA_NOEXCEPT { return 0; }
	virtual int numOutletPorts() const CASEMA_NOEXCEPT { return 1; }
	virtual bool hasInlet() const CASEMA_NOEXCEPT { return false; }
	virtual bool hasOutlet() const CASEMA_NOEXCEPT { return true; }
	virtual bool canAccumulate() const CASEMA_NOEXCEPT { return true; }

	static const char* identifier() { return "INLET"; }
	virtual const char* unitOperationName() const CASEMA_NOEXCEPT { return "INLET"; }

	virtual bool configure(io::IParameterProvider& paramProvider);
	virtual void setFlowRates(mpfr::mpreal const* in, mpfr::mpreal const* out) { }
	virtual void setSectionTimes(mpfr::mpreal const* secTimes, int nTimes) CASEMA_NOEXCEPT;

	virtual void evaluate(const mpfr::mpcomplex& s, Eigen::Ref<Eigen::MatrixCmp> h, Eigen::Ref<Eigen::VectorCmp> g) const CASEMA_NOEXCEPT;
	virtual void evaluate(const mpfr::mpcomplex& s, const mpfr::mpreal& z, Eigen::Ref<Eigen::MatrixCmp> h, Eigen::Ref<Eigen::VectorCmp> g) const CASEMA_NOEXCEPT;

	virtual bool hasValidEstimate() const CASEMA_NOEXCEPT;
	virtual mpfr::mpreal estimate(const mpfr::mpreal& abscissa) const CASEMA_NOEXCEPT;
	virtual mpfr::mpreal timeDomainUpperBound(mpfr::mpreal const* in) const CASEMA_NOEXCEPT;
	virtual mpfr::mpreal truncationError(const mpfr::mpreal& M, const mpfr::mpreal& abscissa, const mpfr::mpreal& T, const mpfr::mpreal& nTerms) const CASEMA_NOEXCEPT { return -1; }
	virtual mpfr::mpreal inverseTruncationError(const mpfr::mpreal& M, const mpfr::mpreal& abscissa, const mpfr::mpreal& T, const mpfr::mpreal& error) const CASEMA_NOEXCEPT { return 0; }

protected:

	std::vector<mpfr::mpreal> _sectionTimes;

	std::vector<mpfr::mpreal> _const;
	std::vector<mpfr::mpreal> _lin;
	std::vector<mpfr::mpreal> _quad;
	std::vector<mpfr::mpreal> _cub;

private:

	template <typename eval_t>
	eval_t evalSection(int i, const eval_t& s) const CASEMA_NOEXCEPT;
};

} // namespace model
} // namespace casema

#endif  // CASEMA_INLETMODEL_HPP_
