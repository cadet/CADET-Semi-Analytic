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
 * Defines the model system, which connects multiple unit operations.
 */

#ifndef CASEMA_MODELSYSTEM_HPP_
#define CASEMA_MODELSYSTEM_HPP_

#include "MPReal.hpp"
#include "MPComplex.hpp"
#include "MPComplexEigen.hpp"
#include "SlicedVector.hpp"

#include <vector>
#include <functional>

namespace casema
{

namespace io
{
	class IParameterProvider;
}

namespace model
{

class UnitOperation;

/**
 * @brief Defines a system of unit operations models
 * @details
 */
class ModelSystem
{
public:

	struct Workspace
	{
		Eigen::MatrixCmp tempMat;
		Eigen::MatrixCmp sysMat;
		Eigen::VectorCmp sysVec;
	};

	ModelSystem();
	~ModelSystem() CASEMA_NOEXCEPT;

	void addModel(UnitOperation* unitOp);
	int numModels() const CASEMA_NOEXCEPT;
	int numOutlets() const CASEMA_NOEXCEPT;
	UnitOperation const* unitOperation(int i) const CASEMA_NOEXCEPT { return _models[i]; }
	UnitOperation* unitOperation(int i) CASEMA_NOEXCEPT { return _models[i]; }

	bool configure(io::IParameterProvider& paramProvider);
	void setSectionTimes(mpfr::mpreal const* secTimes, int nTimes) CASEMA_NOEXCEPT;

	void evaluate(const mpfr::mpcomplex& s, mpfr::mpcomplex* res, Workspace& ws) const CASEMA_NOEXCEPT;

	bool hasValidEstimate() const CASEMA_NOEXCEPT;
	mpfr::mpreal timeDomainUpperBound() const CASEMA_NOEXCEPT;
	mpfr::mpreal truncationError(const mpfr::mpreal& abscissa, const mpfr::mpreal& T, const mpfr::mpreal& nTerms) const CASEMA_NOEXCEPT;
	mpfr::mpreal inverseTruncationError(const mpfr::mpreal& abscissa, const mpfr::mpreal& T, const mpfr::mpreal& error) const CASEMA_NOEXCEPT;
	std::vector<mpfr::mpreal> inverseTruncationErrorOfUnits(const mpfr::mpreal& abscissa, const mpfr::mpreal& T, const mpfr::mpreal& error) const CASEMA_NOEXCEPT;

	void setBesselZeros(mpfr::mpreal const* zeros, int n) CASEMA_NOEXCEPT;
	bool needsBesselZeros() const CASEMA_NOEXCEPT;

	Workspace makeWorkspace() const;

protected:

	void calcUnitFlowRateCoefficients();
	void assembleMatrix();
	int totalNumInletPorts() const CASEMA_NOEXCEPT;
	int totalNumOutletPorts() const CASEMA_NOEXCEPT;
	void addDefaultPortsToConnectionList(std::vector<double>& conList) const;
	void addDefaultDynamicFlowRatesToConnectionList(std::vector<double>& conList) const;
	void setModelFlowRates();
	void checkConnectionList(const std::vector<double>& conn, std::vector<int>& connOnly, std::vector<double>& flow, std::vector<double>& flowLin, std::vector<double>& flowQuad, std::vector<double>& flowCub, int idxSwitch) const;
	void assembleOffsets();

	mpfr::mpreal evaluateConstantOfPath(int idxNode, const mpfr::mpreal& abscissa, const mpfr::mpreal& ToverPi, const std::function<void(int, const mpfr::mpreal&)>& funcUnit) const CASEMA_NOEXCEPT;

	std::vector<UnitOperation*> _models; //!< Unit operation models
	util::SlicedVector<int> _connections; //!< Vector of connection lists for each section
	util::SlicedVector<mpfr::mpreal> _flowRates; //!< Vector of connection flow rates for each section
	util::SlicedVector<mpfr::mpreal> _flowRatesLin; //!< Vector of linear coefficients of connection flow rates for each section
	util::SlicedVector<mpfr::mpreal> _flowRatesQuad; //!< Vector of quadratic coefficients of connection flow rates for each section
	util::SlicedVector<mpfr::mpreal> _flowRatesCub; //!< Vector of cubic coefficients of connection flow rates for each section
	std::vector<int> _switchSectionIndex; //!< Holds indices of sections where valves are switched
	util::SlicedVector<int> _linearModelOrdering; //!< Dependency-consistent ordering of unit operation models for linear execution (for each switch)
	bool _hasDynamicFlowRates;

	util::SlicedVector<int> _reverseAdjList;
	std::vector<int> _terminalIdx;
	std::vector<int> _promotedTerminalIdx;

	util::SlicedVector<mpfr::mpreal> _totalInletFlow; //!< Total flow rate into each inlet at the current section
	util::SlicedVector<mpfr::mpreal> _totalOutletFlow; //!< Total flow rate into each outlet at the current section

	std::vector<int> _offsetUnitInlet;
	std::vector<int> _offsetUnitOutlet;
	Eigen::MatrixCmp _flowMatQt;
	Eigen::MatrixXmp _flowMatQtR;
};

} // namespace model
} // namespace casema

#endif  // CASEMA_MODELSYSTEM_HPP_
