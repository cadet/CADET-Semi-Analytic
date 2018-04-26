// =============================================================================
//  CADET-semi-analytic - The semi analytic extension of
//  		CADET - The Chromatography Analysis and Design Toolkit
//  
//  Copyright © 2015-2018: Samuel Leweke¹
//                                      
//    ¹ Forschungszentrum Juelich GmbH, IBG-1, Juelich, Germany.
//  
//  All rights reserved. This program and the accompanying materials
//  are made available under the terms of the GNU Public License v3.0 (or, at
//  your option, any later version) which accompanies this distribution, and
//  is available at http://www.gnu.org/licenses/gpl.html
// =============================================================================

#ifndef IOEXCEPTION_HPP_
#define IOEXCEPTION_HPP_

#include <string>
#include <stdexcept>

namespace casema
{

class IOException : public std::runtime_error
{
public:
	IOException(const std::string& message) : std::runtime_error(message) { }
};

}  // namespace casema


#endif /* IOEXCEPTION_HPP_ */
