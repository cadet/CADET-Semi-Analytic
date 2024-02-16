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

#include "model/ModelSystem.hpp"
#include "model/UnitOperation.hpp"
#include "io/ParameterProvider.hpp"
#include "Exceptions.hpp"
#include "GraphAlgos.hpp"

//#include <Eigen/QR>
#include <Eigen/LU>

#include <sstream>
#include <iomanip>
#include <cstring>
#include <algorithm>

namespace
{
	template <class Elem_t>
	inline bool contains(const typename std::vector<Elem_t>& vec, const Elem_t& item)
	{
		const typename std::vector<Elem_t>::const_iterator it = std::find(vec.begin(), vec.end(), item);
		return it != vec.end();
	}

	/**
	 * @brief Find array index of unit operation with given id
	 * @details Unit operation index does not need to match index of unit operation in an array.
	 * @param [in] models List of models
	 * @param [in] unitOpIdx Unit operation index to look for
	 * @return Array index of the unit operation identified by unit operation index or invalid array index if the unit operation has not been found
	 */
	inline int indexOfUnitOp(const std::vector<casema::model::UnitOperation*>& models, int unitOpIdx)
	{
		for (int i = 0; i < models.size(); ++i)
		{
			if (models[i]->unitOperationId() == unitOpIdx)
				return i;
		}
		return models.size();
	}

	/**
	 * @brief Returns whether a given unit operation is a terminal node in the network
	 * @param [in] conn Connection list
	 * @param [in] size Number of rows in the connection list
	 * @param [in] unitOpIdx Index of the unit operation to check
	 * @return @c true if the unit operation is a terminal node in the network, @c false otherwise
	 */
	inline bool isTerminal(int const* const conn, int size, int unitOpIdx)
	{
		for (int i = 0; i < size; ++i)
		{
			if (conn[6*i] == unitOpIdx)
				return false;
		}
		return true;
	}

	/**
	 * @brief Returns whether one or more unit operations have multiple ports
	 * @param [in] models List of unit operation models
	 * @return @c true if the list contains at least one unit operation model with multiple ports, @c false otherwise
	 */
	inline bool hasMultiPortUnits(const std::vector<casema::model::UnitOperation*>& models)
	{
		for (casema::model::UnitOperation const* m : models)
		{
			if ((m->numInletPorts() > 1) || (m->numOutletPorts() > 1))
				return true;
		}
		return false;
	}

	template <typename T>
	std::string toSciString(const T val, const int prec = 6)
	{
		std::ostringstream out;
		out << std::scientific << std::setprecision(prec) << val;
		return out.str();
	}
}

namespace casema
{

namespace model
{

ModelSystem::ModelSystem()
{
}

ModelSystem::~ModelSystem() CASEMA_NOEXCEPT
{
	for (UnitOperation* model : _models)
		delete model;
}

void ModelSystem::addModel(UnitOperation* unitOp)
{
	// Check for unique unit operation id
	if (indexOfUnitOp(_models, unitOp->unitOperationId()) < _models.size())
		throw InvalidParameterException("Cannot add model because of already existing unit operation id " + std::to_string(unitOp->unitOperationId()));

	_models.push_back(unitOp);
}

int ModelSystem::numModels() const CASEMA_NOEXCEPT
{
	return _models.size();
}

int ModelSystem::numOutlets() const CASEMA_NOEXCEPT
{
	return _offsetUnitOutlet.back();
}

int ModelSystem::totalNumInletPorts() const CASEMA_NOEXCEPT
{
	int nPorts = 0;
	for (UnitOperation const* m : _models)
		nPorts += m->numInletPorts();

	return nPorts;
}

int ModelSystem::totalNumOutletPorts() const CASEMA_NOEXCEPT
{
	int nPorts = 0;
	for (UnitOperation const* m : _models)
		nPorts += std::max(m->numOutletPorts(), 1);

	return nPorts;
}

bool ModelSystem::configure(io::IParameterProvider& paramProvider)
{
	bool forceSystemSolution = false;
	if (paramProvider.exists("solver"))
	{
		paramProvider.pushScope("solver");
		if (paramProvider.exists("LINEAR_SOLUTION_MODE"))
			forceSystemSolution = paramProvider.getInt("LINEAR_SOLUTION_MODE") == 1;
		paramProvider.popScope();
	}

	// Read connections of unit operations
	paramProvider.pushScope("connections");

	const int numSwitches = paramProvider.getInt("NSWITCHES");
	
	// TODO Improve those very conservative size estimates
	_switchSectionIndex.clear();
	_switchSectionIndex.reserve(numSwitches);
	_connections.clear();
	_connections.reserve(numSwitches * 6 * _models.size() * _models.size(), numSwitches);
	_flowRates.clear();
	_flowRates.reserve(numSwitches * _models.size() * _models.size(), numSwitches);
	_linearModelOrdering.reserve(numSwitches * _models.size(), numSwitches);
	_linearModelOrdering.clear();

	// Default: ports are not given in connection list
	bool conListHasPorts = false;

	// Override default by user option
	if (paramProvider.exists("CONNECTIONS_INCLUDE_PORTS"))
		conListHasPorts = paramProvider.getBool("CONNECTIONS_INCLUDE_PORTS");

	// If we have unit operations with multiple ports, we require ports in connection list
	if (hasMultiPortUnits(_models))
		conListHasPorts = true;

	// Default: no dynamic flow rates
	_hasDynamicFlowRates = false;

	// Override default by user option
	if (paramProvider.exists("CONNECTIONS_INCLUDE_DYNAMIC_FLOW"))
		_hasDynamicFlowRates = paramProvider.getBool("CONNECTIONS_INCLUDE_DYNAMIC_FLOW");

	std::ostringstream oss;
	for (int i = 0; i < numSwitches; ++i)
	{
		oss.str("");
		oss << "switch_" << std::setfill('0') << std::setw(3) << std::setprecision(0) << i;

		paramProvider.pushScope(oss.str());

		_switchSectionIndex.push_back(paramProvider.getInt("SECTION"));

		if ((i > 0) && (_switchSectionIndex.back() <= _switchSectionIndex[_switchSectionIndex.size() - 2]))
			throw InvalidParameterException("SECTION index has to be monotonically increasing (switch " + std::to_string(i) + ")");

		std::vector<double> connFlow = paramProvider.getDoubleArray("CONNECTIONS");
		if (!conListHasPorts)
		{
			if (_hasDynamicFlowRates)
			{
				if ((connFlow.size() % 8) != 0)
					throw InvalidParameterException("CONNECTIONS matrix has to have 8 columns if CONNECTIONS_INCLUDE_PORTS is disabled and CONNECTIONS_INCLUDE_DYNAMIC_FLOW is enabled");
			}
			else
			{
				if ((connFlow.size() % 5) != 0)
					throw InvalidParameterException("CONNECTIONS matrix has to have 5 columns if CONNECTIONS_INCLUDE_PORTS is disabled and CONNECTIONS_INCLUDE_DYNAMIC_FLOW is disabled");
			}

			addDefaultPortsToConnectionList(connFlow);
		}

		if (!_hasDynamicFlowRates)
		{
			if ((connFlow.size() % 7) != 0)
				throw InvalidParameterException("CONNECTIONS matrix has to have 5 or 7 columns if CONNECTIONS_INCLUDE_DYNAMIC_FLOW is disabled");

			addDefaultDynamicFlowRatesToConnectionList(connFlow);
		}

		if ((connFlow.size() % 10) != 0)
			throw InvalidParameterException("CONNECTIONS matrix has to have 10 columns");

		std::vector<int> conn(connFlow.size() / 10 * 6, 0);
		std::vector<double> fr(connFlow.size() / 10, 0.0);
		std::vector<double> frLin(connFlow.size() / 10, 0.0);
		std::vector<double> frQuad(connFlow.size() / 10, 0.0);
		std::vector<double> frCub(connFlow.size() / 10, 0.0);

		checkConnectionList(connFlow, conn, fr, frLin, frQuad, frCub, i);

		_connections.pushBackSlice(conn);

		// Convert double to mpreal while pushing into the SlicedVector
		// also register parameter to enable sensitivities
		if (fr.size() > 0)
		{
			_flowRates.pushBack(fr[0]);
			_flowRatesLin.pushBack(frLin[0]);
			_flowRatesQuad.pushBack(frQuad[0]);
			_flowRatesCub.pushBack(frCub[0]);
			for (int j = 1; j < fr.size(); ++j)
			{
				_flowRates.pushBackInLastSlice(fr[j]);
				_flowRatesLin.pushBackInLastSlice(frLin[j]);
				_flowRatesQuad.pushBackInLastSlice(frQuad[j]);
				_flowRatesCub.pushBackInLastSlice(frCub[j]);

				// Check if a previous identical connection (except for component indices) exists
				bool found = false;
				for (int k = 0; k < j; ++k)
				{
					if ((conn[6*k] == conn[6*j]) && (conn[6*k+1] == conn[6*j+1]) && (conn[6*k+2] == conn[6*j+2]) && (conn[6*k+3] == conn[6*j+3]))
					{
						found = true;
						break;
					}
				}
			}
		}
		else
		{
			// Add empty slice
			_flowRates.pushBackSlice(nullptr, 0);
			_flowRatesLin.pushBackSlice(nullptr, 0);
			_flowRatesQuad.pushBackSlice(nullptr, 0);
			_flowRatesCub.pushBackSlice(nullptr, 0);
		}

		// Auto detect solution method
		const util::SlicedVector<int> adjList = graph::adjacencyListFromConnectionList(conn.data(), _models.size(), conn.size() / 6);
		std::vector<int> topoOrder(0);
		const bool hasCycles = graph::topologicalSort(adjList, topoOrder);

		if (hasCycles || forceSystemSolution)
			_linearModelOrdering.pushBackSlice(0);
		else
			_linearModelOrdering.pushBackSlice(topoOrder);

		paramProvider.popScope();
	}

	paramProvider.popScope();

	if (_switchSectionIndex[0] != 0)
		throw InvalidParameterException("First element of SECTION in connections group has to be 0");

	if (numSwitches > 1)
		throw InvalidParameterException("Only one switch is supported");

	if (_hasDynamicFlowRates)
	{
		for (int i = 0; i < _flowRatesLin.size(); ++i)
		{
			if ((_flowRatesLin.native(i) != 0.0) || (_flowRatesQuad.native(i) != 0.0) || (_flowRatesCub.native(i) != 0.0))
				throw InvalidParameterException("Dynamic flow rates unsupported");
		}
	}

	// Apply flow rates
	calcUnitFlowRateCoefficients();
	setModelFlowRates();
	assembleMatrix();

	// Prepare reverse adjacency list for backwards traversal
	_reverseAdjList = std::move(graph::reverseAdjacencyListFromConnectionList(_connections[0], _models.size(), _connections.sliceSize(0) / 6));

	// Find indices of terminal units in the network
	const util::SlicedVector<int> adjList = graph::adjacencyListFromConnectionList(_connections[0], _models.size(), _connections.sliceSize(0) / 6);
	for (int i = 0; i < _models.size(); ++i)
	{
		if (adjList.sliceSize(i) == 0)
		{
			_terminalIdx.push_back(i);
			if (std::strcmp(_models[i]->unitOperationName(), "OUTLET") == 0)
			{
				// Promote predecessors to terminal nodes
				int const* const idxPre = _reverseAdjList[i];
				for (int j = 0; j < _reverseAdjList.sliceSize(i); ++j)
				{
					if (!contains(_promotedTerminalIdx, idxPre[j]))
						_promotedTerminalIdx.push_back(idxPre[j]);
				}
			}
			else
				_promotedTerminalIdx.push_back(i);
		}
	}

	return true;
}

/**
 * @brief Add default ports to connection list
 * @details Adds source and destination ports of @c -1 to the connection list. The list is
 *          overwritten, that is, all pointers and iterators are invalidated.
 * @param [in,out] conList On entry, list of connections without ports;
 *                         on exit, list of connections including default ports
 */
void ModelSystem::addDefaultPortsToConnectionList(std::vector<double>& conList) const
{
	const int strideSrc = _hasDynamicFlowRates ? 8 : 5;
	const int strideDest = _hasDynamicFlowRates ? 10 : 7;
	const int numRows = conList.size() / strideSrc;
	std::vector<double> newList(numRows * strideDest, -1.0);

	double const* src = conList.data();
	double* dest = newList.data();
	for (int i = 0; i < numRows; ++i, src += strideSrc, dest += strideDest)
	{
		// Always copy unit operation indices
		dest[0] = src[0];
		dest[1] = src[1];

		// Skip ports (default filled with -1)

		// Copy the rest of the line
		std::copy_n(src + 2, strideDest - 4, dest + 4);
	}

	conList = std::move(newList);
}


/**
 * @brief Add default dynamic flow rates to connection list
 * @details Adds linear, quadratic, and cubic flow rate coefficients to the list
 *          The list is overwritten, that is, all pointers and iterators are invalidated.
 * @param [in,out] conList On entry, list of connections without dynamic flow rates;
 *                         on exit, list of connections including dynamic flow rates
 */
void ModelSystem::addDefaultDynamicFlowRatesToConnectionList(std::vector<double>& conList) const
{
	const int numRows = conList.size() / 7;
	std::vector<double> newList(numRows * 10, 0.0);

	double const* src = conList.data();
	double* dest = newList.data();
	for (int i = 0; i < numRows; ++i, src += 7, dest += 10)
	{
		std::copy_n(src, 7, dest);
	}

	conList = std::move(newList);
}

/**
 * @brief Checks the given unit operation connection list and reformats it
 * @details Throws an exception if something is incorrect. Reformats the connection list by
 *          substituting unit operation IDs with local indices.
 * @param [in] conn Matrix with 7 columns holding all connections. The matrix is expected
 *             to be in row-major storage format.
 * @param [out] connOnly Matrix with 6 columns holding all connections. While @p conn contains
 *              connection indices and flow rate, this matrix only holds connection indices.
 *              It is assumed to be pre-allocated (same number of rows as @p conn). The unit
 *              operation IDs are substituted by the corresponding indices of the unit operations
 *              in the local _models vector.
 * @param [out] flowRates Vector with flow rates for each connection in the list. It is assumed
 *              to be pre-allocated (same number of rows as @p conn).
 * @param [out] flowRatesLin Vector with linear flow rate coefficients for each connection in the list. It is assumed
 *              to be pre-allocated (same number of rows as @p conn).
 * @param [out] flowRatesQuad Vector with quadratic flow rate coefficients for each connection in the list. It is assumed
 *              to be pre-allocated (same number of rows as @p conn).
 * @param [out] flowRatesCub Vector with cubic flow rate coefficients for each connection in the list. It is assumed
 *              to be pre-allocated (same number of rows as @p conn).
 * @param [in] idxSwitch Index of the valve switch that corresponds to this connection list
 */
void ModelSystem::checkConnectionList(const std::vector<double>& conn, std::vector<int>& connOnly, std::vector<double>& flowRates, std::vector<double>& flowRatesLin, std::vector<double>& flowRatesQuad, std::vector<double>& flowRatesCub, int idxSwitch) const
{
	std::vector<double> totalInflow(_models.size(), 0.0);
	std::vector<double> totalInflowLin(_models.size(), 0.0);
	std::vector<double> totalInflowQuad(_models.size(), 0.0);
	std::vector<double> totalInflowCub(_models.size(), 0.0);
	std::vector<double> totalOutflow(_models.size(), 0.0);
	std::vector<double> totalOutflowLin(_models.size(), 0.0);
	std::vector<double> totalOutflowQuad(_models.size(), 0.0);
	std::vector<double> totalOutflowCub(_models.size(), 0.0);

	for (int i = 0; i < conn.size() / 10; ++i)
	{
		// Extract current connection
		int uoSource = static_cast<int>(conn[10*i]);
		int uoDest = static_cast<int>(conn[10*i+1]);
		const int portSource = static_cast<int>(conn[10*i+2]);
		const int portDest = static_cast<int>(conn[10*i+3]);
		const int compSource = static_cast<int>(conn[10*i+4]);
		const int compDest = static_cast<int>(conn[10*i+5]);
		double fr = conn[10*i + 6];
		double frLin = conn[10*i + 7];
		double frQuad = conn[10*i + 8];
		double frCub = conn[10*i + 9];

		if (uoSource < 0)
			throw InvalidParameterException("In CONNECTIONS matrix (switch " + std::to_string(idxSwitch) + " row " + std::to_string(i) + "): Source unit operation id has to be at least 0 in connection");
		if (uoDest < 0)
			throw InvalidParameterException("In CONNECTIONS matrix (switch " + std::to_string(idxSwitch) + " row " + std::to_string(i) + "): Destination unit operation id has to be at least 0 in connection");

		// Convert to index
		uoSource = indexOfUnitOp(_models, uoSource);
		uoDest = indexOfUnitOp(_models, uoDest);

		if (static_cast<int>(uoSource) >= _models.size())
			throw InvalidParameterException("In CONNECTIONS matrix (switch " + std::to_string(idxSwitch) + " row " + std::to_string(i) + "): Source unit operation id " + std::to_string(uoSource) + " not found in connection");
		if (static_cast<int>(uoDest) >= _models.size())
			throw InvalidParameterException("In CONNECTIONS matrix (switch " + std::to_string(idxSwitch) + " row " + std::to_string(i) + "): Destination unit operation id " + std::to_string(uoDest) + " not found in connection");

		// Check if unit operations have inlets and outlets
		if (!_models[uoSource]->hasOutlet())
			throw InvalidParameterException("In CONNECTIONS matrix (switch " + std::to_string(idxSwitch) + " row " + std::to_string(i) + "): Source unit operation " + std::to_string(_models[uoSource]->unitOperationId()) + " does not have an outlet");
		if (!_models[uoDest]->hasInlet())
			throw InvalidParameterException("In CONNECTIONS matrix (switch " + std::to_string(idxSwitch) + " row " + std::to_string(i) + "): Source unit operation " + std::to_string(_models[uoDest]->unitOperationId()) + " does not have an inlet");

		// Check port indices
		if ((portSource >= 0) && (static_cast<int>(portSource) >= _models[uoSource]->numOutletPorts()))
			throw InvalidParameterException("In CONNECTIONS matrix (switch " + std::to_string(idxSwitch) + " row " + std::to_string(i) + "): Source port index (" + std::to_string(portSource) + ") exceeds number of outlet ports (" + std::to_string(_models[uoSource]->numOutletPorts()) + ")");
		if ((portDest >= 0) && (static_cast<int>(portDest) >= _models[uoDest]->numInletPorts()))
			throw InvalidParameterException("In CONNECTIONS matrix (switch " + std::to_string(idxSwitch) + " row " + std::to_string(i) + "): Destination port index (" + std::to_string(portDest) + ") exceeds number of inlet ports (" + std::to_string(_models[uoDest]->numInletPorts()) + ")");

		if (((portSource < 0) && (portDest >= 0)) || ((portSource >= 0) && (portDest < 0)))
			throw InvalidParameterException("In CONNECTIONS matrix (switch " + std::to_string(idxSwitch) + " row " + std::to_string(i) + "): Only source or destination (not both) are set to connect all ports in connection from unit " 
				+ std::to_string(_models[uoSource]->unitOperationId()) + " to " + std::to_string(_models[uoDest]->unitOperationId()));

		if ((portSource < 0) && (portDest < 0) && (_models[uoSource]->numOutletPorts() != _models[uoDest]->numInletPorts()))
			throw InvalidParameterException("In CONNECTIONS matrix (switch " + std::to_string(idxSwitch) + " row " + std::to_string(i) + "): Number of ports not equal when connecting all ports from unit " 
				+ std::to_string(_models[uoSource]->unitOperationId()) + " (" + std::to_string(_models[uoSource]->numOutletPorts()) + " ports) to " + std::to_string(_models[uoDest]->unitOperationId()) + " (" + std::to_string(_models[uoDest]->numInletPorts()) + " ports)");

		// Check component indices
		if ((compSource >= 0) && (static_cast<int>(compSource) >= _models[uoSource]->numComponents()))
			throw InvalidParameterException("In CONNECTIONS matrix (switch " + std::to_string(idxSwitch) + " row " + std::to_string(i) + "): Source component index (" + std::to_string(compSource) + ") exceeds number of components (" + std::to_string(_models[uoSource]->numComponents()) + ")");
		if ((compDest >= 0) && (static_cast<int>(compDest) >= _models[uoDest]->numComponents()))
			throw InvalidParameterException("In CONNECTIONS matrix (switch " + std::to_string(idxSwitch) + " row " + std::to_string(i) + "): Destination component index (" + std::to_string(compDest) + ") exceeds number of components (" + std::to_string(_models[uoDest]->numComponents()) + ")");

		if (((compSource < 0) && (compDest >= 0)) || ((compSource >= 0) && (compDest < 0)))
			throw InvalidParameterException("In CONNECTIONS matrix (switch " + std::to_string(idxSwitch) + " row " + std::to_string(i) + "): Only source or destination (not both) are set to connect all components in connection from unit " 
				+ std::to_string(_models[uoSource]->unitOperationId()) + " to " + std::to_string(_models[uoDest]->unitOperationId()));

		if ((compSource < 0) && (compDest < 0) && (_models[uoSource]->numComponents() != _models[uoDest]->numComponents()))
			throw InvalidParameterException("In CONNECTIONS matrix (switch " + std::to_string(idxSwitch) + " row " + std::to_string(i) + "): Number of components not equal when connecting all components from unit " 
				+ std::to_string(_models[uoSource]->unitOperationId()) + " to " + std::to_string(_models[uoDest]->unitOperationId()));
	
		// Add connection to index matrix
		connOnly[6*i] = uoSource;
		connOnly[6*i+1] = uoDest;
		connOnly[6*i+2] = portSource;
		connOnly[6*i+3] = portDest;
		connOnly[6*i+4] = compSource;
		connOnly[6*i+5] = compDest;

		// Add flow rate of connection to balance

		// Check if such a connection has occurred before (for a different component)
		bool found = false;
		for (int j = 0; j < i; ++j)
		{
			if ((conn[10*j] == uoSource) && (uoDest == conn[10*j+1]) && (conn[10*j+2] == portSource) && (portDest == conn[10*j+3]))
			{
				// Take flow rate that appears first
				fr = conn[10*j+6];
				frLin = conn[10*j+7];
				frQuad = conn[10*j+8];
				frCub = conn[10*j+9];
				found = true;
				break;
			}
		}

		// Total flow rates: Only add flow rate once (not for each component)
		if (!found)
		{
			// Add flow rates to balance
			totalInflow[uoDest] += fr;
			totalOutflow[uoSource] += fr;
			totalInflowLin[uoDest] += frLin;
			totalOutflowLin[uoSource] += frLin;
			totalInflowQuad[uoDest] += frQuad;
			totalOutflowQuad[uoSource] += frQuad;
			totalInflowCub[uoDest] += frCub;
			totalOutflowCub[uoSource] += frCub;
		}

		// Add flow rate to list
		flowRates[i] = fr;
		flowRatesLin[i] = frLin;
		flowRatesQuad[i] = frQuad;
		flowRatesCub[i] = frCub;
	}

	// Check flow rate balance
	for (int i = 0; i < _models.size(); ++i)
	{
		// Unit operations with only one port (inlet or outlet) do not need to balance their flows
		if ((totalInflow[i] >= 0.0) && (totalOutflow[i] == 0.0) && _models[i]->hasInlet() && !_models[i]->hasOutlet())
			continue;
		if ((totalInflow[i] == 0.0) && (totalOutflow[i] >= 0.0) && !_models[i]->hasInlet() && _models[i]->hasOutlet())
			continue;

		// Terminal unit operations do not need to balance their flows
		if ((totalOutflow[i] >= 0.0) && isTerminal(connOnly.data(), connOnly.size() / 6, i))
			continue;

		// Unit operations that can accumulate cannot be checked
		if (_models[i]->canAccumulate())
			continue;

		// Check balance
		const double diff = std::abs(totalInflow[i] - totalOutflow[i]);
		if ((diff >= 1e-15) || (diff > 1e-15 * std::abs(totalOutflow[i])))
			throw InvalidParameterException("In CONNECTIONS matrix (switch " + std::to_string(idxSwitch) + "): Flow rate balance is not closed for unit operation " + std::to_string(i) + ", imbalanced by " + toSciString(diff));

		const double diffLin = std::abs(totalInflowLin[i] - totalOutflowLin[i]);
		if ((diffLin >= 1e-15) || (diffLin > 1e-15 * std::abs(totalOutflowLin[i])))
			throw InvalidParameterException("In CONNECTIONS matrix (switch " + std::to_string(idxSwitch) + "): Linear flow rate coefficient balance is not closed for unit operation " + std::to_string(i) + ", imbalanced by " + toSciString(diffLin));

		const double diffQuad = std::abs(totalInflowQuad[i] - totalOutflowQuad[i]);
		if ((diffQuad >= 1e-15) || (diffQuad > 1e-15 * std::abs(totalOutflowQuad[i])))
			throw InvalidParameterException("In CONNECTIONS matrix (switch " + std::to_string(idxSwitch) + "): Quadratic flow rate coefficient balance is not closed for unit operation " + std::to_string(i) + ", imbalanced by " + toSciString(diffQuad));

		const double diffCub = std::abs(totalInflowCub[i] - totalOutflowCub[i]);
		if ((diffCub >= 1e-15) || (diffCub > 1e-15 * std::abs(totalOutflowCub[i])))
			throw InvalidParameterException("In CONNECTIONS matrix (switch " + std::to_string(idxSwitch) + "): Cubic flow rate coefficient balance is not closed for unit operation " + std::to_string(i) + ", imbalanced by " + toSciString(diffCub));
	}

	// TODO: Check for conflicting entries
	// TODO: Plausibility check of total connections
}

/**
* @brief Calculate inlet and outlet flow rate coefficients for each unit operation in current section
*/
void ModelSystem::calcUnitFlowRateCoefficients()
{
	const int curSwitchIndex = 0;

	// Calculate total flow rate for each inlet
	int const* const ptrConn = _connections[curSwitchIndex];
	mpfr::mpreal const* const ptrRate = _flowRates[curSwitchIndex];

	// Reset total flows back to zero
	_totalInletFlow.reserve(totalNumInletPorts(), _models.size());
	_totalOutletFlow.reserve(totalNumOutletPorts(), _models.size());
	for (UnitOperation const* m : _models)
	{
		_totalInletFlow.pushBackSlice(m->numInletPorts());
		_totalOutletFlow.pushBackSlice(m->numOutletPorts());
	}
	_totalInletFlow.fill(0.0);
	_totalOutletFlow.fill(0.0);

	// Compute total volumetric inflow for each unit operation port
	for (int i = 0; i < _connections.sliceSize(curSwitchIndex) / 6; ++i)
	{
		// Extract current connection
		const int uoSource = ptrConn[6*i];
		const int uoDest = ptrConn[6*i + 1];
		const int portSource = ptrConn[6*i + 2];
		const int portDest = ptrConn[6*i + 3];

		// Check if the same connection has appeared before (with different components)
		bool skip = false;
		for (int j = 0; j < i; ++j)
		{
			if ((ptrConn[6*j] == uoSource) && (ptrConn[6*j + 1] == uoDest) && (ptrConn[6*j + 2] == portSource) && (ptrConn[6*j + 3] == portDest))
			{
				skip = true;
				break;
			}
		}

		// Skip this row in connection list if there was an identical previous connection (except for component indices)
		if (skip)
			continue;

		// Use the first flow rate from uoSource to uoDest

		if (portDest < 0)
		{
			for (int j = 0; j < _models[uoDest]->numInletPorts(); ++j)
			{
				_totalInletFlow(uoDest, j) += ptrRate[i];
			}
		}
		else
		{
			_totalInletFlow(uoDest, portDest) += ptrRate[i];
		}

		if (portSource < 0)
		{
			for (int j = 0; j < _models[uoSource]->numOutletPorts(); ++j)
			{
				_totalOutletFlow(uoSource, j) += ptrRate[i];
			}
		}
		else
		{
			_totalOutletFlow(uoSource, portSource) += ptrRate[i];
		}
	}
}

/**
 * @brief Updates inlet and outlet flow rates of the given unit operation
 * @details Updates the corresponding slice of _flowRateIn and _flowRateOut.
 * @param[in] t Time
 * @param[in] idxUnit Unit operation index
 */
void ModelSystem::setModelFlowRates()
{
	for (int i = 0; i < _models.size(); ++i)
	{
		mpfr::mpreal* const in = _totalInletFlow[i];
		mpfr::mpreal* const out = _totalOutletFlow[i];

		_models[i]->setFlowRates(in, out);
	}
}

void ModelSystem::assembleOffsets()
{
	_offsetUnitInlet.clear();
	_offsetUnitOutlet.clear();
	_offsetUnitInlet.reserve(_models.size());
	_offsetUnitOutlet.reserve(_models.size());

	int maxInlet = 0;
	int maxOutlet = 0;
	int curOffsetInlet = 0;
	int curOffsetOutlet = 0;
	for (UnitOperation const* uo : _models)
	{
		_offsetUnitInlet.push_back(curOffsetInlet);
		curOffsetInlet += uo->numInletPorts() * uo->numComponents();

		_offsetUnitOutlet.push_back(curOffsetOutlet);
		curOffsetOutlet += std::max(uo->numOutletPorts(), 1) * uo->numComponents();
	}

	_offsetUnitInlet.push_back(curOffsetInlet);
	_offsetUnitOutlet.push_back(curOffsetOutlet);
}

ModelSystem::Workspace ModelSystem::makeWorkspace() const
{
	int maxInlet = 0;
	int maxOutlet = 0;
	for (UnitOperation const* uo : _models)
	{
		maxInlet = std::max(maxInlet, uo->numInletPorts() * uo->numComponents());
		maxOutlet = std::max(maxOutlet, std::max(uo->numOutletPorts(), 1) * uo->numComponents());
	}

	ModelSystem::Workspace ws;
	ws.tempMat = Eigen::MatrixCmp(maxOutlet, maxInlet);
	ws.sysMat = Eigen::MatrixCmp(totalNumOutletPorts(), totalNumOutletPorts());
	ws.sysVec = Eigen::VectorCmp(totalNumOutletPorts());

	return ws;
}


/**
 * @brief Assembles the bottom macro row handling the connections
 * @details Computes flow rates and ratios for coupling unit operations.
 * @param[in] t Simulation time
 */
void ModelSystem::assembleMatrix()
{
	assembleOffsets();

	_flowMatQt = Eigen::MatrixCmp::Zero(totalNumInletPorts(), totalNumOutletPorts());
	_flowMatQtR = Eigen::MatrixXmp::Zero(totalNumInletPorts(), totalNumOutletPorts());

	const double switchStartTime = 0.0;
	const int curSwitchIndex = 0;

	int const* const ptrConn = _connections[curSwitchIndex];
	mpfr::mpreal const* const ptrRate = _flowRates[curSwitchIndex];

	// Bottom macro-row
	// FN
	for (int i = 0; i < _connections.sliceSize(curSwitchIndex) / 6; ++i)
	{
		// Extract current connection
		const int uoSource = ptrConn[6*i];
		const int uoDest = ptrConn[6*i + 1];
		const int portSource = ptrConn[6*i + 2];
		const int portDest = ptrConn[6*i + 3];
		const int compSource = ptrConn[6*i + 4];
		const int compDest = ptrConn[6*i + 5];

		// Obtain index of first connection from uoSource to uoDest
		int idx = i;
		for (int j = 0; j < i; ++j)
		{
			if ((ptrConn[6*j] == uoSource) && (ptrConn[6*j + 1] == uoDest) && (ptrConn[6*j + 2] == portSource) && (ptrConn[6*j + 3] == portDest))
			{
				idx = j;
				break;
			}
		}

		// idx contains the index of the first connection from uoSource to uoDest
		// Hence, ptrRate[idx] is the flow rate to use for this connection

		UnitOperation const* const modelSource = _models[uoSource];
		UnitOperation const* const modelDest = _models[uoDest];

		// The outlet column is the outlet index + component number * outlet stride

		if (portSource == -1)
		{
			for (int j = 0; j < modelSource->numOutletPorts(); ++j)
			{
				const mpfr::mpreal totInFlow = _totalInletFlow(uoDest, j);

				// Ignore ports with incoming flow rate 0
				if (totInFlow <= 0.0)
					continue;
				
				const mpfr::mpreal inFlow = ptrRate[idx] / totInFlow;

				if (compSource == -1)
				{

					// Connect all components with the same flow rate
					for (int comp = 0; comp < modelSource->numComponents(); ++comp)
					{
						const int row = _offsetUnitInlet[uoDest] + j * modelDest->numComponents() + comp;
						const int col = _offsetUnitOutlet[uoSource] + j * modelSource->numComponents() + comp;
						_flowMatQt(row, col) = inFlow;
						_flowMatQtR(row, col) = inFlow;
					}
				}
				else
				{
					const int row = _offsetUnitInlet[uoDest] + j * modelDest->numComponents() + compDest;
					const int col = _offsetUnitOutlet[uoSource] + j * modelSource->numComponents() + compSource;
					_flowMatQt(row, col) = inFlow;
					_flowMatQtR(row, col) = inFlow;
				}
			}
		}
		else
		{
			const mpfr::mpreal totInFlow = _totalInletFlow(uoDest, portDest);

			// Ignore ports with incoming flow rate 0
			if (totInFlow <= 0.0)
				continue;

			const mpfr::mpreal inFlow = ptrRate[idx] / totInFlow;

			if (compSource == -1)
			{
				// Connect all components with the same flow rate
				for (int comp = 0; comp < modelSource->numComponents(); ++comp)
				{
					const int row = _offsetUnitInlet[uoDest] + portDest * modelDest->numComponents() + comp;
					const int col = _offsetUnitOutlet[uoSource] + portSource * modelSource->numComponents() + comp;
					_flowMatQt(row, col) = inFlow;
					_flowMatQtR(row, col) = inFlow;
				}
			}
			else
			{
				const int row = _offsetUnitInlet[uoDest] + portDest * modelDest->numComponents() + compDest;
				const int col = _offsetUnitOutlet[uoSource] + portSource * modelSource->numComponents() + compSource;
				_flowMatQt(row, col) = inFlow;
				_flowMatQtR(row, col) = inFlow;
			}
		}
	}
}

void ModelSystem::evaluate(const mpfr::mpcomplex& s, mpfr::mpcomplex* res, Workspace& ws) const CASEMA_NOEXCEPT
{
	if (_linearModelOrdering.sliceSize(0) == 0)
	{
		// Network is cyclic

		// Apply block diagonal matrix of multipliers
		for (int i = 0; i < _models.size(); ++i)
		{
			UnitOperation const* const m = _models[i];
			const int numIn = m->numInletPorts() * m->numComponents();
			const int numOut = std::max(m->numOutletPorts(), 1) * m->numComponents();

			m->evaluate(s, ws.tempMat.block(0, 0, numOut, numIn), ws.sysVec.middleRows(_offsetUnitOutlet[i], numOut));
			ws.sysMat.middleRows(_offsetUnitOutlet[i], numOut) = ws.tempMat.block(0, 0, numOut, numIn) * _flowMatQt.middleRows(_offsetUnitInlet[i], numIn);
		}

		ws.sysMat -= Eigen::MatrixCmp::Identity(ws.sysMat.rows(), ws.sysMat.cols());

		Eigen::Map<Eigen::VectorCmp> result(res, ws.sysMat.cols());
		result = ws.sysMat.fullPivLu().solve(-ws.sysVec);
	}
	else
	{
		for (int i = 0; i < _offsetUnitOutlet.back(); ++i)
			res[i] = 0.0;

		Eigen::Map<Eigen::VectorCmp> result(res, _offsetUnitOutlet.back());

		// Network is acyclic
		int const* const order = _linearModelOrdering[0];
		for (int i = _linearModelOrdering.sliceSize(0) - 1; i >= 0; --i)
		{
			const int curIdx = order[i];
			UnitOperation const* const m = _models[curIdx];
			const int numIn = m->numInletPorts() * m->numComponents();
			const int numOut = std::max(m->numOutletPorts(), 1) * m->numComponents();

			m->evaluate(s, ws.tempMat.block(0, 0, numOut, numIn), ws.sysVec.middleRows(_offsetUnitOutlet[curIdx], numOut));

			if (numIn == 0)
				result.middleRows(_offsetUnitOutlet[curIdx], numOut) = ws.sysVec.middleRows(_offsetUnitOutlet[curIdx], numOut);
			else
				result.middleRows(_offsetUnitOutlet[curIdx], numOut) = ws.tempMat.block(0, 0, numOut, numIn) * _flowMatQt.middleRows(_offsetUnitInlet[curIdx], numIn) * result + ws.sysVec.middleRows(_offsetUnitOutlet[curIdx], numOut);
		}
	}
}

bool ModelSystem::hasValidEstimate() const CASEMA_NOEXCEPT
{
	// Check if network is a DAG
	if (_linearModelOrdering.sliceSize(0) == 0)
		return false;

	for (UnitOperation const* unit : _models)
	{
		if (!unit->hasValidEstimate())
			return false;
	}

	// Check if inlet is directly connected to outlet
	for (int i : _promotedTerminalIdx)
	{
		if (std::strcmp(_models[i]->unitOperationName(), "INLET") == 0)
			return false;
	}

	return true;
}

mpfr::mpreal ModelSystem::timeDomainUpperBound() const CASEMA_NOEXCEPT
{
	mpfr::mpreal v(0.0);
	Eigen::VectorXmp inBound = Eigen::VectorXmp::Zero(totalNumInletPorts());
	Eigen::VectorXmp outBound = Eigen::VectorXmp::Zero(totalNumOutletPorts());

	int const* const order = _linearModelOrdering[0];
	for (int i = _linearModelOrdering.sliceSize(0) - 1; i >= 0; --i)
	{
		const int curIdx = order[i];
		UnitOperation const* const m = _models[curIdx];
		const int numIn = m->numInletPorts() * m->numComponents();
		const int numOut = m->numOutletPorts() * m->numComponents();

		if (numIn > 0)
			inBound.middleRows(_offsetUnitInlet[curIdx], numIn) = _flowMatQtR.middleRows(_offsetUnitInlet[curIdx], numIn) * outBound;

		const mpfr::mpreal ob = m->timeDomainUpperBound(inBound.data() + _offsetUnitInlet[curIdx]);
		outBound.middleRows(_offsetUnitOutlet[curIdx], numOut).array() = ob;

		v = max(v, ob);
	}

	return v;
}

mpfr::mpreal ModelSystem::truncationError(const mpfr::mpreal& abscissa, const mpfr::mpreal& T, const mpfr::mpreal& nTerms) const CASEMA_NOEXCEPT
{
	std::vector<mpfr::mpreal> truncErrOfNode(_models.size(), 0);
	const mpfr::mpreal ToverPi = T / mpfr::const_pi();
	for (int i : _promotedTerminalIdx)
	{
		int const* const prev = _reverseAdjList[i];
		const int s = _reverseAdjList.sliceSize(i);

		mpfr::mpreal M(0);
		for (int j = 0; j < s; ++j)
			M += evaluateConstantOfPath(prev[j], abscissa, ToverPi, [&](int idxNode, const mpfr::mpreal& nodeM)
			{
				truncErrOfNode[idxNode] = max(truncErrOfNode[idxNode], _models[idxNode]->truncationError(nodeM, abscissa, T, nTerms));
			}) * _flowMatQtR(_offsetUnitInlet[i], _offsetUnitOutlet[prev[j]]);

		truncErrOfNode[i] = max(truncErrOfNode[i], _models[i]->truncationError(M, abscissa, T, nTerms));
	}

	return *std::max_element(truncErrOfNode.begin(), truncErrOfNode.end());
}

mpfr::mpreal ModelSystem::inverseTruncationError(const mpfr::mpreal& abscissa, const mpfr::mpreal& T, const mpfr::mpreal& error) const CASEMA_NOEXCEPT
{
	const std::vector<mpfr::mpreal> nTermsOfNode = inverseTruncationErrorOfUnits(abscissa, T, error);
	return *std::max_element(nTermsOfNode.begin(), nTermsOfNode.end());
}

std::vector<mpfr::mpreal> ModelSystem::inverseTruncationErrorOfUnits(const mpfr::mpreal& abscissa, const mpfr::mpreal& T, const mpfr::mpreal& error) const CASEMA_NOEXCEPT
{
	const mpfr::mpreal ToverPi = T / mpfr::const_pi();
	std::vector<mpfr::mpreal> nTermsOfNode(_models.size(), 0);
	for (int i : _promotedTerminalIdx)
	{
		int const* const prev = _reverseAdjList[i];
		const int s = _reverseAdjList.sliceSize(i);

		mpfr::mpreal M(0);
		for (int j = 0; j < s; ++j)
			M += evaluateConstantOfPath(prev[j], abscissa, ToverPi, [&](int idxNode, const mpfr::mpreal& nodeM)
			{
				nTermsOfNode[idxNode] = max(nTermsOfNode[idxNode], ceil(_models[idxNode]->inverseTruncationError(nodeM, abscissa, T, error)));
			}) * _flowMatQtR(_offsetUnitInlet[i], _offsetUnitOutlet[prev[j]]);

		nTermsOfNode[i] = max(nTermsOfNode[i], ceil(_models[i]->inverseTruncationError(M, abscissa, T, error)));
	}

	return nTermsOfNode;
}

mpfr::mpreal ModelSystem::evaluateConstantOfPath(int idxNode, const mpfr::mpreal& abscissa, const mpfr::mpreal& ToverPi, const std::function<void(int, const mpfr::mpreal&)>& funcUnit) const CASEMA_NOEXCEPT
{
	const int s = _reverseAdjList.sliceSize(idxNode);
	if (s == 0)
	{
		// This is a terminal node (e.g., inlet)
		funcUnit(idxNode, 0);
		return _models[idxNode]->estimate(abscissa) * ToverPi;
	}

	mpfr::mpreal out(0);
	int const* const prev = _reverseAdjList[idxNode];
	for (int i = 0; i < s; ++i)
		out += evaluateConstantOfPath(prev[i], abscissa, ToverPi, funcUnit) * _flowMatQtR(_offsetUnitInlet[idxNode], _offsetUnitOutlet[prev[i]]); // evaluateConstantOfPath(prev[i], abscissa) * _flowMatQtR(_offsetUnitInlet[idxNode], _offsetUnitOutlet[prev[i]]);

	funcUnit(idxNode, out);
	return out * _models[idxNode]->estimate(abscissa) * ToverPi;
}


void ModelSystem::setSectionTimes(mpfr::mpreal const* secTimes, int nTimes) CASEMA_NOEXCEPT
{
	for (UnitOperation* m : _models)
		m->setSectionTimes(secTimes, nTimes);
}

void ModelSystem::setBesselZeros(mpfr::mpreal const* zeros, int n) CASEMA_NOEXCEPT
{
	for (UnitOperation* m : _models)
		m->setBesselZeros(zeros, n);
}

bool ModelSystem::needsBesselZeros() const CASEMA_NOEXCEPT
{
	for (UnitOperation const* m : _models)
	{
		if (m->needsBesselZeros())
			return true;
	}
	return false;
}

}  // namespace model

}  // namespace casema
