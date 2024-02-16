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
 * Defines exceptions.
 */

#ifndef CASEMA_EXCEPTIONS_HPP_
#define CASEMA_EXCEPTIONS_HPP_

#include <stdexcept>

namespace casema
{

/**
 * @brief Signals invalid parameter or option values
 */
class InvalidParameterException : public std::domain_error
{
public:
	explicit InvalidParameterException(const std::string& what_arg) : std::domain_error(what_arg) { }
	explicit InvalidParameterException(const char* what_arg) : std::domain_error(what_arg) { }
};


} // namespace casema

#endif  // CASEMA_EXCEPTIONS_HPP_
