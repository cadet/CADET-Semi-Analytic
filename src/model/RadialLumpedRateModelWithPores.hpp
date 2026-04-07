// =============================================================================
//  CADET-semi-analytic - The semi-analytic extension of CADET
//
//  Copyright © 2015-present: Samuel Leweke¹² and the CADET-Semi-Analytic
//  authors, see the AUTHORS file
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
 * Defines the radial flow lumped rate model with pores (rLRMP).
 */

#ifndef CASEMA_RADIALLUMPEDRATEMODELWITHPORES_HPP_
#define CASEMA_RADIALLUMPEDRATEMODELWITHPORES_HPP_

#include "model/LumpedRateModelWithPores.hpp"

namespace casema
{

namespace model
{

/**
 * @brief Radial flow lumped rate model of liquid column chromatography with pores
 * @details See @cite Gu1991, @cite Besselink2013
 *
 * Uses a volume-equivalent coordinate transformation to map the radial
 * flow PDE onto the axial LRMP in the Laplace domain.
 */
class RadialLumpedRateModelWithPores : public LumpedRateModelWithPores
{
public:

	RadialLumpedRateModelWithPores(int unitOpIdx);
	virtual ~RadialLumpedRateModelWithPores() CASEMA_NOEXCEPT;

	static const char* identifier() { return "RADIAL_LUMPED_RATE_MODEL_WITH_PORES"; }
	virtual const char* unitOperationName() const CASEMA_NOEXCEPT { return "RADIAL_LUMPED_RATE_MODEL_WITH_PORES"; }

	virtual bool configure(io::IParameterProvider& paramProvider);

protected:

	virtual void setVelocityFromFlowRate(mpfr::mpreal const* in);

	mpfr::mpreal _innerRadius; //!< Inner column radius \f$ \rho_c \f$
	mpfr::mpreal _outerRadius; //!< Outer column radius \f$ \rho \f$
	mpfr::mpreal _physColLength; //!< Physical column length \f$ L \f$ (axial)
	int _flowDirection; //!< Flow direction: +1 = outward, -1 = inward
};

} // namespace model
} // namespace casema

#endif  // CASEMA_RADIALLUMPEDRATEMODELWITHPORES_HPP_