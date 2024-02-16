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

#ifndef LIBCASEMA_IOEXCEPTION_HPP_
#define LIBCASEMA_IOEXCEPTION_HPP_

#include <string>
#include <stdexcept>

namespace casema 
{

namespace io
{

class IOException : public std::runtime_error
{
public:
	IOException(const std::string& message) : std::runtime_error(message) { }
};


} // namespace io

} // namespace casema


#endif /* LIBCASEMA_IOEXCEPTION_HPP_ */
