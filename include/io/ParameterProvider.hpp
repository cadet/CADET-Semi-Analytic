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

#ifndef CASEMA_PARAMPROVIDER_HPP_
#define CASEMA_PARAMPROVIDER_HPP_

#include <string>
#include <vector>
#include <cstdint>

#include "casemaCompilerInfo.hpp"

namespace casema
{
namespace io
{

/**
 * @brief ParameterProvider interface is used for querying parameters and data
 */
class IParameterProvider
{
public:

	virtual ~IParameterProvider() CASEMA_NOEXCEPT { }

	/**
	 * @brief Returns the value of a parameter of type double
	 * 
	 * @param [in] paramName Name of the parameter
	 * @return Parameter value
	 */
	virtual double getDouble(const std::string& paramName) = 0;

	/**
	 * @brief Returns the value of a parameter of type int
	 * 
	 * @param [in] paramName Name of the parameter
	 * @return Parameter value
	 */
	virtual int getInt(const std::string& paramName) = 0;

	/**
	 * @brief Returns the value of a parameter of type uint64_t
	 * 
	 * @param [in] paramName Name of the parameter
	 * @return Parameter value
	 */
	virtual uint64_t getUint64(const std::string& paramName) = 0;

	/**
	 * @brief Returns the value of a parameter of type bool
	 * 
	 * @param [in] paramName Name of the parameter
	 * @return Parameter value
	 */
	virtual bool getBool(const std::string& paramName) = 0;

	/**
	 * @brief Returns the value of a parameter of type string
	 * 
	 * @param [in] paramName Name of the parameter
	 * @return Parameter value
	 */
	virtual std::string getString(const std::string& paramName) = 0;

	/**
	 * @brief Returns a parameter array of type double
	 * 
	 * @param [in] paramName Name of the parameter
	 * @return Parameter values
	 */
	virtual std::vector<double> getDoubleArray(const std::string& paramName) = 0;

	/**
	 * @brief Returns a parameter array of type int
	 * 
	 * @param [in] paramName Name of the parameter
	 * @return Parameter values
	 */
	virtual std::vector<int> getIntArray(const std::string& paramName) = 0;

	/**
	 * @brief Returns a parameter array of type uint64_t
	 * 
	 * @param [in] paramName Name of the parameter
	 * @return Parameter values
	 */
	virtual std::vector<uint64_t> getUint64Array(const std::string& paramName) = 0;

	/**
	 * @brief Returns a parameter array of type bool
	 * 
	 * @param [in] paramName Name of the parameter
	 * @return Parameter values
	 */
	virtual std::vector<bool> getBoolArray(const std::string& paramName) = 0;

	/**
	 * @brief Returns a parameter array of type string
	 * 
	 * @param [in] paramName Name of the parameter
	 * @return Parameter values
	 */
	virtual std::vector<std::string> getStringArray(const std::string& paramName) = 0;

	/**
	 * @brief Checks whether a given parameter exists
	 * 
	 * @param [in] paramName Name of the parameter
	 * @return @c true if the parameter exists, otherwise @c false
	 */
	virtual bool exists(const std::string& paramName) = 0;

	/**
	 * @brief Checks whether a given parameter is an array
	 * 
	 * @param [in] paramName Name of the parameter
	 * @return @c true if the parameter is an array, otherwise @c false
	 */
	virtual bool isArray(const std::string& paramName) = 0;

	/**
	 * @brief Returns the number of elements (of an array) in a given field
	 * @param [in] paramName Name of the parameter
	 * @return Number of elements in the given field
	 */
	virtual std::size_t numElements(const std::string& paramName) = 0;

	/**
	 * @brief Changes to a given namespace subscope
	 *
	 * @param [in] scope Name of the scope
	 */
	virtual void pushScope(const std::string& scope) = 0;

	/**
	 * @brief Changes to the parent of the current namespace scope
	 */
	virtual void popScope() = 0;
};

} // namespace io
} // namespace casema

#endif  // CASEMA_PARAMPROVIDER_HPP_
