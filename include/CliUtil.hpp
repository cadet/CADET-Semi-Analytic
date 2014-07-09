// =============================================================================
//  CADET-semi-analytic - The semi analytic extension of
//  		CADET - The Chromatography Analysis and Design Toolkit
//  
//  Copyright © 2015: Samuel Leweke¹
//                                      
//    ¹ Forschungszentrum Juelich GmbH, IBG-1, Juelich, Germany.
//  
//  All rights reserved. This program and the accompanying materials
//  are made available under the terms of the GNU Public License v3.0 (or, at
//  your option, any later version) which accompanies this distribution, and
//  is available at http://www.gnu.org/licenses/gpl.html
// =============================================================================

#ifndef CASEMA_CLIUTIL_HPP_
#define CASEMA_CLIUTIL_HPP_

#include <string>
#include <mpreal.h>

namespace TCLAP
{
    // Enable TCLAP support of mpreal datatypes
    template<>
    struct ArgTraits<mpfr::mpreal>
    {
        typedef StringLike ValueCategory;
    };
    
    template<>
    void SetString(mpfr::mpreal &dst, const std::string &src)
    {
        const std::size_t prec = mpfr::mpreal::get_default_prec();
        const std::size_t inprec = mpfr::digits2bits(src.size()+2);
        
        mpfr::mpreal x(0, std::max(prec, inprec));
        mpfr_strtofr (x.mpfr_ptr(), src.c_str(), 0, 10, mpfr::mpreal::get_default_rnd());

        dst = x;
    }
}

#endif
