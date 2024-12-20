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
 * Provides tools for managing scopes of IParameterProvider
 */

#ifndef CASEMA_PARAMREADERSCOPES_HPP_
#define CASEMA_PARAMREADERSCOPES_HPP_

#include "io/ParameterProvider.hpp"
#include "Exceptions.hpp"

#include <string>
#include <sstream>
#include <iomanip>

namespace casema
{
	/**
	 * @brief ParameterProvider scope guard for a multiplexed group
	 * @details Opens and closes the ParameterProvider scope.
	 */
	class MultiplexedScopeSelector
	{
	public:

		/**
		 * @brief Enters a multiplexed scope in the IParameterProvider
		 * @details If multiplexing is disabled, the scope @p baseName or @c baseName_000 is entered, if it exists.
		 *          If the scope does not exist, but it is @p required, an exception is thrown.
		 *          If multiplxing is enabled, the scope @c baseName_xxx is entered, where @c xxx is a given @p index.
		 *          If this scope does not exist and @p fallbackSingle is enabled, we additionally check the scope @p baseName.
		 *          If all this fails, an exception is thrown if the scope is @p required.
		 *          
		 *          On destruction of the object, the scope is popped back from the IParameterProvider (if it has been entered).
		 * @param [in,out] paramProvider IParameterProvider instance
		 * @param [in] baseName Base name of the group
		 * @param [in] single Determines whether a single scope is multiplexed on multiple items
		 * @param [in] index Index of the non-multiplexed scope
		 * @param [in] fallbackSingle Determines whether a fallback scope (@p baseName) is tried if the scope identified by @p index is not found
		 * @param [in] required Determines whether the scope is required
		 */
		MultiplexedScopeSelector(io::IParameterProvider& paramProvider, const std::string& baseName, bool single, unsigned int index, bool fallbackSingle, bool required) : _pp(paramProvider)
		{
			_active = true;
			if (single)
			{
				if (paramProvider.exists(baseName))
					paramProvider.pushScope(baseName);
				else if (paramProvider.exists(baseName + "_000"))
					paramProvider.pushScope(baseName + "_000");
				else
				{
					_active = false;
					if (required)
						throw InvalidParameterException("Group \"" + baseName + "\" or \"" + baseName + "_000\" required");
				}
			}
			else
			{
				std::ostringstream oss;
				oss << baseName + "_" << std::setfill('0') << std::setw(3) << std::setprecision(0) << index;

				if (paramProvider.exists(oss.str()))
					paramProvider.pushScope(oss.str());
				else
				{
		 			// If there is just one type, allow legacy "baseName" scope
					if (paramProvider.exists(baseName) && fallbackSingle)
						paramProvider.pushScope(baseName);
					else
					{
						_active = false;
						if (required)
							throw InvalidParameterException("Group \"" + oss.str() + "\" or \"" + baseName + "\" required");
					}
				}
			}
		}

		/**
		 * @brief Enters a multiplexed scope in the IParameterProvider
		 * @details The scope @p baseName or @c baseName_000 is entered, if it exists.
		 *          If the scope does not exist, but it is @p required, an exception is thrown.
		 *          
		 *          On destruction of the object, the scope is popped back from the IParameterProvider (if it has been entered).
		 * @param [in,out] paramProvider IParameterProvider instance
		 * @param [in] baseName Base name of the group
		 * @param [in] required Determines whether the scope is required
		 */
		MultiplexedScopeSelector(io::IParameterProvider& paramProvider, const std::string& baseName, bool required) : _pp(paramProvider)
		{
			_active = true;
			if (paramProvider.exists(baseName))
				paramProvider.pushScope(baseName);
			else if (paramProvider.exists(baseName + "_000"))
				paramProvider.pushScope(baseName + "_000");
			else
			{
				_active = false;
				if (required)
					throw InvalidParameterException("Group \"" + baseName + "\" or \"" + baseName + "_000\" required");
			}
		}

		/**
		 * @brief Enters a multiplexed scope in the IParameterProvider
		 * @details The scope @c baseName_xxx is entered, where @c xxx is a given @p index.
		 *          If this scope does not exist and @p fallbackSingle is enabled, we additionally check the scope @p baseName.
		 *          If all this fails, an exception is thrown if the scope is @p required.
		 *          
		 *          On destruction of the object, the scope is popped back from the IParameterProvider (if it has been entered).
		 * @param [in,out] paramProvider IParameterProvider instance
		 * @param [in] baseName Base name of the group
		 * @param [in] index Index of the non-multiplexed scope
		 * @param [in] fallbackSingle Determines whether a fallback scope (@p baseName) is tried if the scope identified by @p index is not found
		 * @param [in] required Determines whether the scope is required
		 */
		MultiplexedScopeSelector(io::IParameterProvider& paramProvider, const std::string& baseName, unsigned int index, bool fallbackSingle, bool required) : _pp(paramProvider)
		{
			_active = true;
			std::ostringstream oss;
			oss << baseName + "_" << std::setfill('0') << std::setw(3) << std::setprecision(0) << index;

			if (paramProvider.exists(oss.str()))
				paramProvider.pushScope(oss.str());
			else
			{
	 			// If there is just one type, allow legacy "baseName" scope
				if (paramProvider.exists(baseName) && fallbackSingle)
					paramProvider.pushScope(baseName);
				else
				{
					_active = false;
					if (required)
						throw InvalidParameterException("Group \"" + oss.str() + "\" or \"" + baseName + "\" required");
				}
			}
		}

		~MultiplexedScopeSelector()
		{
			if (_active)
			{
				_pp.popScope();
			}
		}

		/**
		 * @brief Returns whether a scope has been entered
		 * @return @c true if a scope has been entered, otherwise @c false
		 */
		inline bool isActive() const CASEMA_NOEXCEPT { return _active; }

	private:
		bool _active;
		io::IParameterProvider& _pp;
	};

} // namespace casema

#endif  // CASEMA_PARAMREADERSCOPES_HPP_
